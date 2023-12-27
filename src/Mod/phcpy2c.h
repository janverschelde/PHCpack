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

static PyObject *py2c_corecount
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of cores available for multithreading. */

/* The wrapping of functions with prototypes in phcpack.h starts here. */

static PyObject *py2c_PHCpack_version_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the version string of PHCpack.
 *   The version string is 40 characters long. */

static PyObject *py2c_set_seed
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Takes the value of the integer given on input
 *   and sets the seed for the random number generators.
 *   This fixing of the seed enables reproducible runs. */

static PyObject *py2c_get_seed
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_define_output_file
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_write_dobldobl_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target system as stored in double double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_quaddobl_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target system as stored in quad double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_standard_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start system as stored in standard double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_dobldobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start system as stored in double double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_quaddobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start system as stored in double double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_read_standard_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file,
 *   in standard double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_standard_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start Laurent system in standard double precision. */

static PyObject *py2c_read_standard_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file,
 *   in standard double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_standard_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target Laurent system in standard double precision. */

static PyObject *py2c_read_dobldobl_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file,
 *   in double double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_dobldobl_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start Laurent system in double double precision. */

static PyObject *py2c_read_dobldobl_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file,
 *   in double double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_dobldobl_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target Laurent system in double double precision. */

static PyObject *py2c_read_quaddobl_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file,
 *   in quad double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_quaddobl_start_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start Laurent system in quad double precision. */

static PyObject *py2c_read_quaddobl_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file,
 *   in quad double precision.
 *   If available on file, also its solutions will be read and stored. */

static PyObject *py2c_write_quaddobl_target_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target Laurent system in quad double precision. */

/* copying systems from and to containers */

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

static PyObject *py2c_copy_standard_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in standard double precision
 *  from the container to the start system. */

static PyObject *py2c_copy_dobldobl_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in double double precision
 *  from the container to the start system. */

static PyObject *py2c_copy_quaddobl_Laurent_container_to_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in quad double precision
 *  from the container to the start system. */

static PyObject *py2c_copy_standard_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in standard double precision
 *  from the container to the target system. */

static PyObject *py2c_copy_dobldobl_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in double double precision
 *  from the container to the target system. */

static PyObject *py2c_copy_quaddobl_Laurent_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Copies the Laurent system in quad double precision
 *  from the container to the target system. */

static PyObject *py2c_copy_standard_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start Laurent system in standard double precision
 *   to the systems container for Laurent systems. */

static PyObject *py2c_copy_dobldobl_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start Laurent system in double double precision
 *   to the systems container for Laurent systems. */

static PyObject *py2c_copy_quaddobl_Laurent_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start Laurent system in quad double precision
 *   to the systems container for Laurent systems. */

static PyObject *py2c_copy_standard_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target Laurent system in standard double precision
 *   to the systems container for Laurent systems. */

static PyObject *py2c_copy_dobldobl_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target Laurent system in double double precision
 *   to the systems container for Laurent systems. */

static PyObject *py2c_copy_quaddobl_Laurent_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target Laurent system in quad double precision
 *   to the systems container for Laurent systems. */

/* creation of homotopy and tracking of all solution paths */

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
 *   On input are two doubles and one positive integer:
 *   (1) the real and imaginary part of the gamma constant;
 *   (2) the power of t in the homotopy.
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
 *   On input are two doubles and one positive integer:
 *   (1) the real and imaginary part of the gamma constant;
 *   (2) the power of t in the homotopy.
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
 *   On input are two doubles and one positive integer:
 *   (1) the real and imaginary part of the gamma constant;
 *   (2) the power of t in the homotopy.
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
 *   On input are two doubles and one positive integer:
 *   (1) the real and imaginary part of the gamma constant;
 *   (2) the power of t in the homotopy.
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

static PyObject *py2c_write_start_solutions
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_get_value_of_continuation_parameter
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current value of a continuation parameter.
 *   On input is an index of a continuation parameter, in the range 1..34,
 *   on return is the current value of that continuation parameter. */

static PyObject *py2c_set_value_of_continuation_parameter
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the value of a continuation parameter.
 *   On input is an index of a continuation parameter, in the range 1..34,
 *   and the new value for the corresponding continuation parameter. */

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

static PyObject *py2c_solve_by_standard_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in standard double precision
 *   to solve a Laurent system stored in the systems container,
 *   starting at the solutions of a stored Laurent start system.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_solve_by_dobldobl_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in double double precision
 *   to solve a Laurent system stored in the systems container,
 *   starting at the solutions of a stored Laurent start system.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_solve_by_quaddobl_Laurent_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in quad double precision
 *   to solve a Laurent system stored in the systems container,
 *   starting at the solutions of a stored Laurent start system.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_clear_standard_operations_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates all data used by solve_by_standard_homotopy_continuation. */

static PyObject *py2c_clear_dobldobl_operations_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates all data used by solve_by_dobldobl_homotopy_continuation. */

static PyObject *py2c_clear_quaddobl_operations_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates all data used by solve_by_quaddobl_homotopy_continuation. */

static PyObject *py2c_clear_standard_Laurent_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data used to solve Laurent systems by homotopy 
 *   continuation in standard double precision. */

static PyObject *py2c_clear_dobldobl_Laurent_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data used to solve Laurent systems by homotopy
 *   continuation in double double precision. */

static PyObject *py2c_clear_quaddobl_Laurent_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data used to solve Laurent systems by homotopy
 *   continuation in quad double precision. */

/* Wrapping of crude path trackers of the jumpstart library starts here. */

static PyObject *py2c_standard_crude_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *   Tracking happens in standard double precision.
 *   On entry is the verbose parameter which is 1 or 0.  
 *   If 1, then the solution vectors are written to screen, otherwise
 *   the crude tracker stays mute.
 *   On return is the failure code, which is zero when all went well.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in standard
 *   double precision have been initialized in the containers. */

static PyObject *py2c_dobldobl_crude_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *   Tracking happens in double double precision.
 *   On entry is the verbose parameter which is 1 or 0.  
 *   If 1, then the solution vectors are written to screen, otherwise
 *   the crude tracker stays mute.
 *   On return is the failure code, which is zero when all went well.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in double
 *   double precision have been initialized in the containers. */

static PyObject *py2c_quaddobl_crude_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *   Tracking happens in quad double precision.
 *   On entry is the verbose parameter which is 1 or 0.  
 *   If 1, then the solution vectors are written to screen, otherwise
 *   the crude tracker stays mute.
 *   On return is the failure code, which is zero when all went well.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in quad
 *   double precision have been initialized in the containers. */

/* The wrapping of copying solutions from and to containers starts here. */

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

static PyObject *py2c_solve_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   systems with coefficients in standard double precision.
 *   Three integers are expected on input: (1) silent or not (0 or 1);
 *   (2) the number of tasks, if 0, then no multitasking; and
 *   (3) the verbose level.
 *   On return, the container for solutions in standard double precision
 *   contains the solutions to the system in the standard systems container. */

static PyObject *py2c_scan_for_symbols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given on input are two arguments: a number and a string.
 *   The string holds the string representation of a polynomial system,
 *   where each polynomial is terminated by a semi colon.
 *   The first argument on input is the number of characters in the string.
 *   On return is the number of symbols used as variables in the system.
 *   This function helps to determine whether a system is square or not. */

static PyObject *py2c_solve_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   systems with coefficients in double double precision.
 *   Three integers are expected on input: (1) silent or not (0 or 1);
 *   (2) the number of tasks, if 0, the no multitasking; and
 *   (3) the verbose level.
 *   On return, the container for solutions in double double precision
 *   contains the solutions to the system in the dobldobl systems container. */

static PyObject *py2c_solve_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   systems with coefficients in quad double precision.
 *   Three integers are expected on input: (1) silent or not (0 or 1);
 *   (2) the number of tasks, if 0, then no multitasking; and
 *   (3) the verbose level.
 *   On return, the container for solutions in quad double precision
 *   contains the solutions to the system in the quaddobl systems container. */

static PyObject *py2c_solve_standard_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   Laurent systems with coefficients in standard double precision.
 *   Three integers are expected on input:
 *   1) a boolean flag silent: if 1, then no intermediate output about
 *   the root counts is printed, if 0, then the solver is verbose;
 *   2) the number of tasks: if 0, then no multitasking is applied,
 *   otherwise as many tasks as the number will run; and
 *   3) the verbose level.
 *   On return, the container for solutions in standard double precision
 *   contains the solutions to the system in the standard Laurent systems
 *   container. */

static PyObject *py2c_solve_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   Laurent systems with coefficients in double double precision.
 *   Two integers are expected on input:
 *   1) a boolean flag silent: if 1, then no intermediate output about
 *   the root counts is printed, if 0, then the solver is verbose;
 *   2) the number of tasks: if 0, then no multitasking is applied,
 *   otherwise as many tasks as the number will run; and
 *   3) the verbose level.
 *   On return, the container for solutions in double double precision
 *   contains the solutions to the system in the dobldobl Laurent systems
 *   container. */

static PyObject *py2c_solve_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   Laurent systems with coefficients in quad double precision.
 *   Two integers are expected on input:
 *   1) a boolean flag silent: if 1, then no intermediate output about
 *   the root counts is printed, if 0, then the solver is verbose;
 *   2) the number of tasks: if 0, then no multitasking is applied,
 *   otherwise as many tasks as the number will run; and
 *   3) the verbose level.
 *   On return, the container for solutions in quad double precision
 *   contains the solutions to the system in the quaddobl Laurent systems
 *   container. */

static PyObject *py2c_set_gamma_constant ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores the gamma constant for later retrieval.
 *   Four parameters are expected on input, two doubles and two integers.
 *   The two doubles are the real and imaginary parts of the gamma.
 *   The two integers are the precision, 1, 2, or 4, respectively for
 *   double, double double, or quad double; and the verbose level. */

static PyObject *py2c_get_gamma_constant ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the gamma constant used by the solve functions.
 *   Two integer parameters are expected on input:
 *   (1) for the precision, 1, 2, or 4, respectively for double,
 *       double double, or quad double precision; and
 *   (2) the verbose level.
 *   The function returns a tuple of two doubles,
 *   for the real and imaginary part of the gamma constant. */

static PyObject *py2c_mixed_volume
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the mixed volume, and the stable mixed volume as well if
 *   the input parameter equals 1.  On return is the mixed volume, or
 *   a tuple with the mixed volume and the stable mixed volume.
 *   The Ada translation of the MixedVol algorithm is applied.
 *   The system in the standard systems container is considered first
 *   and if that container is empty, then the system in the standard
 *   Laurent systems container is taken as input.
 *   A regular mixed-cell configuration is in the cells container. */

static PyObject *py2c_mixed_volume_by_demics 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls DEMiCs to compute the mixed volume of the system in the
 *   standard systems container.  If the standard systems container
 *   is empty, then the system in the standard Laurent systems
 *   container is taken as input.  Returns the mixed volume.
 *   A regular mixed-cell configuration is in the cells container.
 *   The above is for the case if the input parameter equals 0.
 *   If the input parameter equals 1, then on return is a tuple,
 *   which contains the mixed volume and the stable mixed volume. */

static PyObject *py2c_standard_deflate
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in standard double precision to the system and
 *   the solutions stored in the containers.
 *   There are five parameters, two integers and three doubles:
 *   1) the maximum number of iterations per root,
 *   2) the maximum number of deflations per root,
 *   3) tolerance on the error of each root,
 *   4) tolerance on the residual of each root,
 *   5) tolerance on the numerical rank of the Jacobians at each root.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_deflate
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in double double precision to the system and
 *   the solutions stored in the containers.
 *   There are five parameters, two integers and three doubles:
 *   1) the maximum number of iterations per root,
 *   2) the maximum number of deflations per root,
 *   3) tolerance on the error of each root,
 *   4) tolerance on the residual of each root,
 *   5) tolerance on the numerical rank of the Jacobians at each root.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_deflate
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in quad double precision to the system and
 *   the solutions stored in the containers.
 *   There are five parameters, two integers and three doubles:
 *   1) the maximum number of iterations per root,
 *   2) the maximum number of deflations per root,
 *   3) tolerance on the error of each root,
 *   4) tolerance on the residual of each root,
 *   5) tolerance on the numerical rank of the Jacobians at each root.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_Newton_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in standard double precision to the system in
 *   the standard systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Newton_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in double double precision to the system in
 *   the dobldobl systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Newton_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in quad double precision to the system in
 *   the quaddobl systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_multprec_Newton_step
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_standard_condition_report
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the system and solutions in the containers in double precision,
 *   computes a condition report.  On input are the following:
 *   1) maximum number of Newton iterations per solution;
 *   2) tolerance on the residual;
 *   3) tolerance on the forward error;
 *   4) tolerance on the inverse condition number for singularities;
 *   5) a string with the name of the output file,
 *   this string may be empty if no output to file is needed;
 *   6) a verbose flag, either 1 or 0.
 *   On return are the counts of number of solutions that are 
 *   regular, singular, real, complex, clustered, or failures;
 *   along with the frequency tables for the forward errors,
 *   residuals and estimates for the inverse condition numbers.  */

/* The wrapping of functions with prototypes in unisolvers.h starts here. */

static PyObject *py2c_usolve_standard
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_usolve_dobldobl
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_usolve_quaddobl
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_usolve_multprec
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_giftwrap_planar
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the giftwrapping algorithm to a planar point configuration.
 *   On input are an integer and a string:
 *   1) the number of points in the list;
 *   2) the string representation of a Python list of tuples.
 *   On return is the string representation of the vertex points,
 *   sorted so that each two consecutive points define an edge. */

static PyObject *py2c_giftwrap_convex_hull
 ( PyObject *self, PyObject *args );
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
 *   the support of the k-th Laurent polynomial in the container.
 *   The index k is given on input as an integer between 1
 *   and the number of Laurent polynomials. */

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

static PyObject *py2c_giftwrap_initial_form
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the Laurent systems container by its initial form.
 *   There are three input parameters:
 *   1) the dimension, number of coordinates in the inner normal;
 *   2) the number of characters in the string representation for the normal;
 *   3) the string representation of the inner normal.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of the functions in syscon.h starts from here */

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

static PyObject *py2c_syscon_random_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Places in the systems container a random polynomial system
 *   with coefficients in standard double precision.
 *   There are five integers as input parameters:
 *   1) n, the number of polynomials and variables;
 *   2) m, the number of monomials per equation;
 *   3) d, the largest degree of each monomial;
 *   4) c, the type of coefficient: 0 if on the complex unit circle,
 *   1, if all coefficients are one, 2, if all coefficients are
 *   random floats in [-1,+1];
 *   5) neq, the number of equations. */

static PyObject *py2c_syscon_dobldobl_random_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Places in the systems container a random polynomial system
 *   with coefficients in double double precision.
 *   There are five integers as input parameters:
 *   1) n, the number of polynomials and variables;
 *   2) m, the number of monomials per equation;
 *   3) d, the largest degree of each monomial;
 *   4) c, the type of coefficient: 0 if on the complex unit circle,
 *   1, if all coefficients are one, 2, if all coefficients are
 *   random floats in [-1,+1];
 *   5) neq, the number of equations. */

static PyObject *py2c_syscon_quaddobl_random_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Places in the systems container a random polynomial system
 *   with coefficients in double double precision.
 *   There are five integers as input parameters:
 *   1) n, the number of polynomials and variables;
 *   2) m, the number of monomials per equation;
 *   3) d, the largest degree of each monomial;
 *   4) c, the type of coefficient: 0 if on the complex unit circle,
 *   1, if all coefficients are one, 2, if all coefficients are
 *   random floats in [-1,+1];
 *   5) neq, the number of equations. */

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

static PyObject *py2c_syscon_retrieve_term
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_syscon_total_degree
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_syscon_standard_drop_variable_by_name
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

static PyObject *py2c_syscon_standard_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the standard double precision container 
 *   with the same Laurent system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syscon_standard_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the standard double precision container 
 *   with the same Laurent system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_dobldobl_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the double double precision container 
 *   with the same Laurent system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syscon_dobldobl_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the double double precision container 
 *   with the same Laurent system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_quaddobl_Laurent_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the quad double precision container 
 *   with the same Laurent system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syscon_quaddobl_Laurent_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system in the quad double precision container 
 *   with the same Laurent system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_standard_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the standard double precision container
 *   with its transformation in 1-homogeneous coordinates.
 *   There is one integer on input.
 *   If 0, then a random linear equation is added,
 *   otherwise, the linear equation z0 - 1 = 0 is added,
 *   where z0 is the extra homogeneous coordinate. */

static PyObject *py2c_syscon_dobldobl_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the double double precision container
 *   with its transformation in 1-homogeneous coordinates.
 *   There is one integer on input.
 *   If 0, then a random linear equation is added,
 *   otherwise, the linear equation z0 - 1 = 0 is added,
 *   where z0 is the extra homogeneous coordinate. */

static PyObject *py2c_syscon_quaddobl_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the quad double precision container
 *   with its transformation in 1-homogeneous coordinates.
 *   There is one integer on input.
 *   If 0, then a random linear equation is added,
 *   otherwise, the linear equation z0 - 1 = 0 is added,
 *   where z0 is the extra homogeneous coordinate. */

static PyObject *py2c_syscon_add_symbol
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Adds a symbol to the table, with name given in the string,
 *   where the number of characters in the name equals the first
 *   integer argument.  The second input parameter is the string.
 *   This symbol represents the last variable added in the homogeneous
 *   coordinate transformation. */

static PyObject *py2c_syscon_standard_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the standard double precision container
 *   by its transformation to affine coordinates, substituting the
 *   value of the last coordinate by one and removing the last equation. */

static PyObject *py2c_syscon_dobldobl_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the double double precision container
 *   by its transformation to affine coordinates, substituting the
 *   value of the last coordinate by one and removing the last equation. */

static PyObject *py2c_syscon_quaddobl_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the quad double precision container
 *   by its transformation to affine coordinates, substituting the
 *   value of the last coordinate by one and removing the last equation. */

/* The wrapping of the functions in tabform.h starts from here */

static PyObject *py2c_tabform_store_standard_tableau
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On input is the tableau form of a polynomial system, given by
 *   1) the number of equations as an integer,
 *   2) the number of equations as an integer,
 *   3) the number of characters in the 4-th string input,
 *   4) the number of terms in each polynomial, given as a string,
 *   the string representation of a list of integers,
 *   5) the number of characters in the 6-th string input,
 *   6) the coefficients of all terms, given as a string,
 *   the string representation of a list of doubles,
 *   each pair of consecutive doubles represents a complex coefficient,
 *   7) the number of characters in the 7-th string input,
 *   8) the exponents of all terms, given as a string,
 *   the string representation of a list of integers,
 *   9) the verbose flag is an integer.
 *   The tableau form is parsed and the container for systems with
 *   standard double precision coefficients is initialized. */

static PyObject *py2c_tabform_load_standard_tableau
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a 5-tuple with the tableau form of the system with
 *   standard double precision coefficients in the container.
 *   On input is the verbose flag as an integer.
 *   The five items in the returned tuple are
 *   1) the number of equations as an integer,
 *   2) the number of equations as an integer,
 *   3) the number of terms in each polynomial, given as a string,
 *   the string representation of a list of integers,
 *   4) the coefficients of all terms, given as a string,
 *   the string representation of a list of doubles,
 *   each pair of consecutive doubles represents a complex coefficient,
 *   5) the exponents of all terms, given as a string,
 *   the string representation of a list of integers. */

/* The wrapping of functions with prototypes in solcon.h starts here. */

static PyObject *py2c_solcon_read_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in standard double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_standard_solutions_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a solution list to be parsed in standard double precision.
 *   Solutions are read from file and stored in the container for
 *   double precision solutions.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_solcon_read_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in double double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_dobldobl_solutions_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a solution list to be parsed in double double precision.
 *   Solutions are read from file and stored in the container for
 *   double double precision solutions.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_solcon_read_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in quad double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_quaddobl_solutions_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a solution list to be parsed in quad double precision.
 *   Solutions are read from file and stored in the container for
 *   quad double precision solutions.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_solcon_read_multprec_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in arbitrary multiprecision.
 *   Returns the failure code, which is zero when all went well. */

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

static PyObject *py2c_solcon_standard_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in standard double precision. */

static PyObject *py2c_solcon_dobldobl_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in double double precision. */

static PyObject *py2c_solcon_quaddobl_one_homogenization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in double double precision. */

static PyObject *py2c_solcon_standard_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in standard double precision. */

static PyObject *py2c_solcon_dobldobl_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in double double precision. */

static PyObject *py2c_solcon_quaddobl_one_affinization
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in quad double precision. */

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

/* The wrapping of the functions with prototypes in celcon.h starts here. */

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

static PyObject *py2c_celcon_is_stable ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns 1 if the stable mixed cells were stored, returns 0 otherwise. */

static PyObject *py2c_celcon_number_of_original_cells
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of original cells, without artificial original. */

static PyObject *py2c_celcon_number_of_stable_cells
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of stable cells. */

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

static PyObject *py2c_celcon_solve_stable_standard_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th stable mixed cell,
 *   using standard double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in standard double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th stable mixed cell. */

static PyObject *py2c_celcon_solve_stable_dobldobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th stable mixed cell,
 *   using double double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in double double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th stable mixed cell. */

static PyObject *py2c_celcon_solve_stable_quaddobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th stable mixed cell,
 *   using quad double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in quad double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th stable mixed cell. */

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

/* The wrapping of functions with prototypes in intcelcon.h follows. */

static PyObject *py2c_intcelcon_read_mixed_cell_configuration
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name, reads a mixed-cell configuration,
 *   which is then stored in the container. */

static PyObject *py2c_intcelcon_write_mixed_cell_configuration
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the mixed-cell configuration in the container to screen. */

static PyObject *py2c_intcelcon_number_of_cells
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of cells in the cell container. */

static PyObject *py2c_intcelcon_type_of_mixture
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the type of mixture
 *   of the support sets.  This string is the string representation
 *   of a Python list of integers. */

static PyObject *py2c_intcelcon_length_of_supports
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of a list of integers where 
 *   each integer contains the number of points in a support. */

static PyObject *py2c_intcelcon_get_lifted_point
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are three integer numbers: n, i, j, where n is the
 *   length of a lifted point, i is the index of a support, and
 *   j is the index of a point.  Note that both i and j start at
 *   one instead of at zero.  Returns the string representation of
 *   the n coordinates of the j-th point in the i-th lifted support. */

static PyObject *py2c_intcelcon_get_inner_normal
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given on input the dimension of the lifted points and the
 *   index of the mixed cell of interest, returns the string
 *   representation of the inner normal of the mixed cell. */

static PyObject *py2c_intcelcon_number_of_points_in_cell
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given is the index to a cell, starting the count at one.
 *   The second integer on input is the number of distinct supports.
 *   On return is the string representation of the number of points
 *   which span each component of the mixed cell. */

static PyObject *py2c_intcelcon_get_point_in_cell 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the n coordinates of 
 *   the k-th point from the j-th list of the i-th cell.
 *   On input are the four integers: n, i, j, k, respectively
 *   the length of the lifted vectors in the supports,
 *   the index to a cell in the container,
 *   the index to a support of the i-th cell, and
 *   the index to a point in the j-th support of the i-th cell. */

static PyObject *py2c_intcelcon_mixed_volume
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the mixed volume of the supports stored
 *   in the cell container. */

static PyObject *py2c_intcelcon_initialize_supports
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the cell container with the number of distinct supports,
 *   this number is given as the one input parameter.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_intcelcon_set_type_of_mixture
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

static PyObject *py2c_intcelcon_append_lifted_point
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a lifted point to the cells container.
 *   There are three input parameters:
 *   1) the dimension of the point;
 *   2) the index of the support to where to append to; and
 *   3) the string representation of the lifted point.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_intcelcon_make_subdivision
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the mixed cells for the lifted points stored in the
 *   container, with respect to the defined type of mixture. */

static PyObject *py2c_intcelcon_clear_mixed_cell_configuration
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the cells container. */

/* The wrapping of functions with prototypes in scalers.h starts from here. */

static PyObject *py2c_scale_standard_system
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_scale_dobldobl_system
 ( PyObject *self, PyObject *args );
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

static PyObject *py2c_scale_quaddobl_system
 ( PyObject *self, PyObject *args );
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

/* The wrapping of functions with prototypes in reducers.h starts here. */

static PyObject *py2c_linear_reduce_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies linear reduction on the coefficient matrix of the system
 *   in the container for standard double precision.
 *   There is one integer parameter: whether to diagonalize or not. */

static PyObject *py2c_linear_reduce_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies linear reduction on the coefficient matrix of the system
 *   in the container for double double precision.
 *   There is one integer parameter: whether to diagonalize or not. */

static PyObject *py2c_linear_reduce_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies linear reduction on the coefficient matrix of the system
 *   in the container for quad double precision.
 *   There is one integer parameter: whether to diagonalize or not. */

static PyObject *py2c_nonlinear_reduce_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies nonlinear reduction on the system in the container 
 *   for standard double precision.
 *   Three integer numbers are expected on input:
 *   (1) the maximum number of equal degree replacements,
 *   (2) the maximum number of computed S-polynomials,
 *   (3) the maximum number of computed R-polynomials.
 *   The system in the standard container is replace by the reduced system.
 *   Three numbers are returned:
 *   (1) the number of equal degree replacements,
 *   (2) the number of computed S-polynomials,
 *   (3) the number of computed R-polynomials. */

/* The wrapping of the functions with prototypes in sweep.h starts here. */

static PyObject *py2c_sweep_define_parameters_numerically
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the indices to the variables that serve as parameters
 *   numerically, that is: via integer indices.
 *   On entry are three integer numbers and a string.
 *   The string is a string representation of a Python list of integers,
 *   The three integers are the number of equations, the number of variables,
 *   and the number of parameters.  The number of variables m includes the
 *   number of parameters.  Then there should be as many as m indices in
 *   the list of integers to define which of the variables are parameters. */

static PyObject *py2c_sweep_define_parameters_symbolically
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the indices to the variables that serve as parameters
 *   symbolically, that is, as names of variables.
 *   For this to work, the symbol table must be initialized.
 *   On entry are four integer numbers and a string.
 *   The four integers are the number of equations, the number of variables,
 *   the number of parameters (the number of variables m includes the
 *   number of parameters), and the number of characters in the string.
 *   The string contains the names of the parameters, separated by one comma.
 *   For this to work, the symbol table must be initialized, e.g.:
 *   via the reading of a polynomial system. */

static PyObject *py2c_sweep_get_number_of_equations
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of equations. */

static PyObject *py2c_sweep_get_number_of_variables
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of variables. */

static PyObject *py2c_sweep_get_number_of_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of parameters. */

static PyObject *py2c_sweep_get_indices_numerically
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the indices of the variables that are parameters,
 *   as the string representation of a Python list of integers. */

static PyObject *py2c_sweep_get_indices_symbolically
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string with the names of the parameters,
 *   each separated by one space. */

static PyObject *py2c_sweep_clear_definitions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the definitions of the parameters. */

static PyObject *py2c_sweep_set_standard_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the start values for the m parameters in standard double precision,
 *   giving on input an integer m and 2*m doubles, with the consecutive 
 *   real and imaginary parts for the start values of all m parameters.
 *   The doubles are given in a string representation of a Python
 *   list of doubles. */

static PyObject *py2c_sweep_set_standard_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the target values for the m parameters in standard double precision,
 *   giving on input an integer m and 2*m doubles, with the consecutive
 *   real and imaginary parts for the target values of all m parameters. */

static PyObject *py2c_sweep_set_dobldobl_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the start values for the m parameters in double double precision,
 *   giving on input an integer m and 4*m doubles, with the consecutive
 *   real and imaginary parts for the start values of all m parameters. */

static PyObject *py2c_sweep_set_dobldobl_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the target values for the m parameters in double double precision,
 *   giving on input an integer m and 4*m doubles, with the consecutive
 *   real and imaginary parts for the target values of all m parameters. */

static PyObject *py2c_sweep_set_quaddobl_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the start values for the m parameters in quad double precision,
 *   giving on input an integer m and 8*m doubles, with the consecutive
 *   real and imaginary parts for the start values of all m parameters. */

static PyObject *py2c_sweep_set_quaddobl_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the target values for the m parameters in quad double precision,
 *   giving on input an integer m and 8*m doubles, with the consecutive
 *   real and imaginary parts for the target values of all m parameters. */

static PyObject *py2c_sweep_get_standard_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in standard double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_get_standard_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Gets the target values for the parameters in standard double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_get_dobldobl_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in double double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_get_dobldobl_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Gets the target values for the parameters in double double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_get_quaddobl_start
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in quad double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_get_quaddobl_target
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the target values for the parameters in quad double precision,
 *   giving on input the number n of doubles that need to be returned.
 *   On return will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters,
 *   stored in the string representation of a Python list of doubles. */

static PyObject *py2c_sweep_standard_complex_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in standard double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in standard double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The first input parameter is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts given in 2 pointers to doubles. */

static PyObject *py2c_sweep_dobldobl_complex_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in double double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in double double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The first input parameter is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts given in 2 pointers to doubles. */

static PyObject *py2c_sweep_quaddobl_complex_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in quad double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in quad double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The first input parameter is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts given in 2 pointers to doubles. */

static PyObject *py2c_sweep_standard_real_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There are no input arguments to this routine.
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

static PyObject *py2c_sweep_dobldobl_real_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There are no input arguments to this routine.
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

static PyObject *py2c_sweep_quaddobl_real_run
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There are no input arguments to this routine.
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

/* The wrapping for the multiplicity structure starts here. */

static PyObject *py2c_standard_multiplicity_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure in standard double precision.
 *   Required is the presence of a polynomial system in the standard
 *   systems container and a solution in the standard solutions container.
 *   The input parameters are two integers and one double:
 *   order : the maximum differentiation order,
 *   verbose : 1 for verbose, 0 for silent, and
 *   tol : tolerance on the numerical rank.
 *   On return is a tuple: the multiplicity and the values
 *   of the Hilbert function. */

static PyObject *py2c_dobldobl_multiplicity_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure in double double precision.
 *   Required is the presence of a polynomial system in the dobldobl
 *   systems container and a solution in the dobldobl solutions container.
 *   The input parameters are two integers and one double:
 *   order : the maximum differentiation order,
 *   verbose : 1 for verbose, 0 for silent, and
 *   tol : tolerance on the numerical rank.
 *   On return is a tuple: the multiplicity and the values
 *   of the Hilbert function. */

static PyObject *py2c_quaddobl_multiplicity_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure in quad double precision.
 *   Required is the presence of a polynomial system in the quaddobl
 *   systems container and a solution in the quaddobl solutions container.
 *   The input parameters are two integers and one double:
 *   order : the maximum differentiation order,
 *   verbose : 1 for verbose, 0 for silent, and
 *   tol : tolerance on the numerical rank.
 *   On return is a tuple: the multiplicity and the values
 *   of the Hilbert function. */

/* The wrapping of the numerical tropisms container starts here. */

static PyObject *py2c_numbtrop_standard_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the numerical tropisms container,
 *   in standard double precision.  The input parameters are
 *   nbt : number of tropisms;
 *   dim : length_of_each tropism;
 *   wnd : winding numbers, as many as nbt;
 *   dir : nbt*dim doubles with the coordinates of the tropisms;
 *   err : errors on the tropisms, as many doubles as the value of nbt.
 *   The numbers in wnd, dir, and err must be given in one string,
 *   as the string representation of a list of doubles.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_dobldobl_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the numerical tropisms container,
 *   in double double precision.  The input parameters are
 *   nbt : number of tropisms;
 *   dim : length_of_each tropism;
 *   wnd : winding numbers, as many as nbt;
 *   dir : 2*nbt*dim doubles with the coordinates of the tropisms;
 *   err : errors on the tropisms, as many doubles as the value of 2*nbt.
 *   The numbers in wnd, dir, and err must be given in one string,
 *   as the string representation of a list of doubles.
 *   On return is the the failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_quaddobl_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the numerical tropisms container,
 *   in quad double precision.  The input parameters are
 *   nbt : number of tropisms;
 *   dim : length_of_each tropism;
 *   wnd : winding numbers, as many as nbt;
 *   dir : 4*nbt*dim doubles with the coordinates of the tropisms;
 *   err : errors on the tropisms, as many doubles as the value of 4*nbt.
 *   The numbers in wnd, dir, and err must be given in one string,
 *   as the string representation of a list of doubles.
 *   On return is the the failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_standard_retrieve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Retrieves all tropisms stored in standard double precision.
 *
 * ON ENTRY :
 *   nbt    number of tropisms;
 *   dim    length_of_each tropism.
 *
 * ON RETURN :
 *   wnd    winding numbers, as many as nbt;
 *   dir    nbt*dim doubles with the coordinates of the tropisms;
 *   err    errors on the tropisms, as many doubles as the value of nbt.
 *   All return parameters are in one string,
 *   the string representation of a list of doubles.
 *   The failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_dobldobl_retrieve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Retrieves all tropisms stored in double double precision.
 *
 * ON ENTRY :
 *   nbt    number of tropisms;
 *   dim    length_of_each tropism.
 *
 * ON RETURN :
 *   wnd    winding numbers, as many as nbt;
 *   dir    2*nbt*dim doubles with the coordinates of the tropisms;
 *   err    errors on the tropisms, as many doubles as the value of 2*nbt.
 *   All return parameters are in one string,
 *   the string representation of a list of doubles.
 *   The failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_quaddobl_retrieve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Retrieves all tropisms stored in quad double precision.
 *
 * ON ENTRY :
 *   nbt    number of tropisms;
 *   dim    length_of_each tropism.
 *
 * ON RETURN :
 *   wnd    winding numbers, as many as nbt;
 *   dir    4*nbt*dim doubles with the coordinates of the tropisms;
 *   err    errors on the tropisms, as many doubles as the value of 4*nbt.
 *   All return parameters are in one string,
 *   the string representation of a list of doubles.
 *   The failure code, which equals zero if all went well. */

static PyObject *py2c_numbtrop_standard_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of tropisms, stored in standard double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_dobldobl_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of tropisms, stored in double double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_quaddobl_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of tropisms, stored in quad double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_standard_dimension
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the dimension of the tropisms, stored in standard double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_dobldobl_dimension
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the dimension of the tropisms, stored in double double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_quaddobl_dimension
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the dimension of the tropisms, stored in quad double
 *   precision, in the numerical tropisms container. */

static PyObject *py2c_numbtrop_store_standard_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores a tropism given in standard double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by standard_size.
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as dim;
 *   err     the error on the tropism.
 *   All dim+1 doubles are given in one string,
 *   the string representation of a list of doubles. */

static PyObject *py2c_numbtrop_store_dobldobl_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores a tropism given in double double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by standard_size.
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as 2*dim;
 *   err     the error on the tropism, two doubles.
 *   All 2*dim+2 doubles are given in one string,
 *   the string representation of a list of doubles. */

static PyObject *py2c_numbtrop_store_quaddobl_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores a tropism given in quad double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by standard_size.
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as 4*dim;
 *   err     the error on the tropism, four doubles.
 *   All 4*dim+4 doubles are given in one string,
 *   the string representation of a list of doubles. */

static PyObject *py2c_numbtrop_standard_retrieve_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns one tropism, stored in standard double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by standard_size.
 *
 * ON RETURN :
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as dim;
 *   err     the error on the tropism.
 *   All dim+1 doubles are returned in one string,
 *   the string representation of a list of doubles. */

static PyObject *py2c_numbtrop_dobldobl_retrieve_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns one tropism, stored in double double precision.
 * 
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by numbtrop_dobldobl_size.
 *
 * ON RETURN :
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as 2*dim;
 *   err     the error on the tropism, two doubles.
 *   All 2*dim+2 doubles are returned in one string,
 *   the string representation of a list of doubles. */

static PyObject *py2c_numbtrop_quaddobl_retrieve_tropism
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns one tropism, stored in quad double precision.
 *
 * ON ENTRY :
 *   dim     the length of the tropism vector;
 *   idx     the index of the tropism, indexing starts at one,
 *           and ends at nbt, what is returned by numbtrop_quaddobl_size.
 *
 * ON RETURN :
 *   wnd     estimated winding number;
 *   dir     coordinates of the tropisms, as many as 4*dim;
 *   err     the error on the tropism, four doubles.
 *   All 4*dim+4 doubles are returned in one string,
 *   the string representation of a list of doubles.*/

static PyObject *py2c_numbtrop_standard_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the stored numerically computed tropisms,
 *   computed in standard double precision. */

static PyObject *py2c_numbtrop_dobldobl_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the stored numerically computed tropisms,
 *   computed in double double precision. */

static PyObject *py2c_numbtrop_quaddobl_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the stored numerically computed tropisms,
 *   computed in quad double precision. */

/* The wrapping of functions with prototypes in witset.h starts here. */

static PyObject *py2c_embed_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the container with its embedding of dimension d.
 *   The dimension d is given as the first integer parameter on input.
 *   The second integer parameter indicates the precision, either 0, 1, or 2,
 *   respectively for double, double double, or quad double precision.
 *   The third integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system with coefficients in standard double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system with coefficients in double double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system with coefficients in quad double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_standard_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system with coefficients in standard double
 *   precision in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system with coefficients in double double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_embed_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the Laurent system with coefficients in quad double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   The second integer parameter is the verbose level.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_swap_symbols_for_standard_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the polynomial system with standard 
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

static PyObject *py2c_swap_symbols_for_dobldobl_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the polynomial system with double
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

static PyObject *py2c_swap_symbols_for_quaddobl_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the polynomial system with quad
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

static PyObject *py2c_swap_symbols_for_standard_Laurent_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the Laurent system with standard 
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

static PyObject *py2c_swap_symbols_for_dobldobl_Laurent_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the Laurent system with double
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

static PyObject *py2c_swap_symbols_for_quaddobl_Laurent_witness_set
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the slack variables in the Laurent system with quad
 *   double precision coefficients and its corresponding solutions in the 
 *   containers so the slack variables appear at the end.  On input are
 *   two integers: the total number of variables; and
 *   the number of slack variables, or the dimension of the set.
 *   This permutation is necessary to consider the system and solutions
 *   stored in containers as a witness set. */

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

static PyObject *py2c_standard_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in standard double precision using the stored
 *   Laurent systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in double double precision using the stored
 *   Laurent systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Laurent_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in quad double precision using the stored
 *   Laurent systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_set_standard_to_mute
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in standard double
 *   precision to silent. */

static PyObject *py2c_factor_set_dobldobl_to_mute
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in double double
 *   precision to silent. */

static PyObject *py2c_factor_set_quaddobl_to_mute
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in quad double
 *   precision to silent. */

static PyObject *py2c_factor_set_standard_to_verbose
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in standard double
 *   precision to verbose. */

static PyObject *py2c_factor_set_dobldobl_to_verbose
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in double double
 *   precision to verbose. */

static PyObject *py2c_factor_set_quaddobl_to_verbose
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations in quad double
 *   precision to verbose. */

static PyObject *py2c_factor_define_output_file_with_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the output file for the factorization.
 *   On input are an integer and a string:
 *   1) the integer equals the number of characters in the string; and
 *   2) the string contains the name of a file.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_standard_assign_labels
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns labels, replacing the multiplicity field of each solution
 *   in standard double precision stored in the container.
 *   On entry are two integers:
 *   1) n, the number of coordinates of the solutions;
 *   2) nbsols, the number of solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_dobldobl_assign_labels
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns labels, replacing the multiplicity field of each solution
 *   in double double precision stored in the container.
 *   On entry are two integers:
 *   1) n, the number of coordinates of the solutions;
 *   2) nbsols, the number of solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_quaddobl_assign_labels
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns labels, replacing the multiplicity field of each solution
 *   in quad double precision stored in the container.
 *   On entry are two integers:
 *   1) n, the number of coordinates of the solutions;
 *   2) nbsols, the number of solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_initialize_standard_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by an ordinary polynomial system in standard double precision.
 *   The embedded system is taken from the polynomial systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_dobldobl_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by an ordinary polynomial system in double double precision.
 *   The embedded system is taken from the polynomial systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_quaddobl_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by an ordinary polynomial system in quad double precision.
 *   The embedded system is taken from the polynomial systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_standard_Laurent_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by a Laurent polynomial system in standard double precision.
 *   The embedded system is taken from the Laurent systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_dobldobl_Laurent_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by a Laurent polynomial system in double double precision.
 *   The embedded system is taken from the Laurent systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_quaddobl_Laurent_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set,
 *   defined by a Laurent polynomial system in quad double precision.
 *   The embedded system is taken from the Laurent systems container
 *   and the generic points are taken from the solutions container.
 *   On entry is the dimension or the number of hyperplanes
 *   to intersect the positive dimensional solution set with. */

static PyObject *py2c_factor_initialize_standard_monodromy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the internal data structures for n loops,
 *   to factor a k-dimensional solution component of degree d,
 *   in standard double precision.
 *   There are three integers on input, in the following order:
 *   1) n, the number of loops;
 *   2) d, the degree of the solution set;
 *   3) k, the dimensional of the solution set.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_factor_initialize_dobldobl_monodromy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the internal data structures for n loops,
 *   to factor a k-dimensional solution component of degree d,
 *   in double double precision.
 *   There are three integers on input, in the following order:
 *   1) n, the number of loops;
 *   2) d, the degree of the solution set;
 *   3) k, the dimensional of the solution set.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_factor_initialize_quaddobl_monodromy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the internal data structures for n loops,
 *   to factor a k-dimensional solution component of degree d,
 *   in quad double precision.
 *   There are three integers on input, in the following order:
 *   1) n, the number of loops;
 *   2) d, the degree of the solution set;
 *   3) k, the dimensional of the solution set.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_factor_standard_trace_grid_diagnostics
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple of two doubles with the diagnostics on the
 *   trace grid computed in standard double precision.
 *   The first double is the largest error of the samples.
 *   The second double is the smallest distance between two samples. */

static PyObject *py2c_factor_dobldobl_trace_grid_diagnostics
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple of two doubles with the diagnostics on the
 *   trace grid computed in double double precision.
 *   The first double is the largest error of the samples.
 *   The second double is the smallest distance between two samples. */

static PyObject *py2c_factor_quaddobl_trace_grid_diagnostics
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple of two doubles with the diagnostics on the
 *   trace grid computed in quad double precision.
 *   The first double is the largest error of the samples.
 *   The second double is the smallest distance between two samples. */

static PyObject *py2c_factor_store_standard_solutions
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Stores the solutions in the container, in standard double precision,
 *   to the data for monodromy loops. */

static PyObject *py2c_factor_store_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Stores the solutions in the container, in double double precision,
 *   to the data for monodromy loops. */

static PyObject *py2c_factor_store_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Stores the solutions in the container, in quad double precision,
 *   to the data for monodromy loops. */

static PyObject *py2c_factor_restore_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Restores the first initialized solutions, in standard double precision,
 *   from sampler to the container. */

static PyObject *py2c_factor_restore_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Restores the first initialized solutions, in double double precision,
 *   from sampler to the container. */

static PyObject *py2c_factor_restore_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Restores the first initialized solutions, in quad double precision,
 *   from sampler to the container. */

static PyObject *py2c_factor_standard_track_paths
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Tracks as many paths as defined by witness set,
 *   in standard double precision.
 *   On input is an integer, which must be 1 if the witness set is
 *   defined by a Laurent polynomial system.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_dobldobl_track_paths
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Tracks as many paths as defined by witness set,
 *   in double double precision.
 *   On input is an integer, which must be 1 if the witness set is
 *   defined by a Laurent polynomial system.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_quaddobl_track_paths
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Tracks as many paths as defined by witness set,
 *   in quad double precision.
 *   On input is an integer, which must be 1 if the witness set is
 *   defined by a Laurent polynomial system.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_swap_standard_slices
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back, in standard double precision.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_swap_dobldobl_slices
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back, in double double precision.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_swap_quaddobl_slices
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back, in quad double precision.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_new_standard_slices
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Generates k random slides in n-space in standard double precision.
 *   The k and the n are the two input parameters.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_new_dobldobl_slices
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Generates k random slides in n-space in double double precision.
 *   The k and the n are the two input parameters.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_new_quaddobl_slices
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Generates k random slides in n-space in quad double precision.
 *   The k and the n are the two input parameters.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_set_standard_trace_slice
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns the constant coefficient of the first slice,
 *   in standard double precision.
 *   On entry is a flag to indicate if it was the first time or not.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_set_dobldobl_trace_slice
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns the constant coefficient of the first slice,
 *   in double double precision.
 *   On entry is a flag to indicate if it was the first time or not.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_set_quaddobl_trace_slice
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns the constant coefficient of the first slice,
 *   in quad double precision.
 *   On entry is a flag to indicate if it was the first time or not.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_store_standard_gammas
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores the gamma constants in standard double precision
 *   for the sampler in the monodromy loops.
 *   Generates as many random complex constants as the value on input.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_store_dobldobl_gammas
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores the gamma constants in double double precision
 *   for the sampler in the monodromy loops.
 *   Generates as many random complex constants as the value on input.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_store_quaddobl_gammas
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores the gamma constants in quad double precision
 *   for the sampler in the monodromy loops.
 *   Generates as many random complex constants as the value on input.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_permutation_after_standard_loop
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For a set of degree d, computes the permutation using the solutions
 *   most recently stored, after a loop in standard double precision.
 *   The number d is the input parameter of this function.
 *   On return is the string representation of the permutation. */

static PyObject *py2c_factor_permutation_after_dobldobl_loop
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For a set of degree d, computes the permutation using the solutions
 *   most recently stored, after a loop in double double precision.
 *   The number d is the input parameter of this function.
 *   On return is the string representation of the permutation. */

static PyObject *py2c_factor_permutation_after_quaddobl_loop
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For a set of degree d, computes the permutation using the solutions
 *   most recently stored, after a loop in quad double precision.
 *   The number d is the input parameter of this function.
 *   On return is the string representation of the permutation. */

static PyObject *py2c_factor_update_standard_decomposition
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Updates the decomposition with the given permutation of d elements,
 *   computed in standard double precision.
 *   On entry are two integers and one string:
 *   1) d, the number of elements in the permutation;
 *   2) nc, the number of characters in the string;
 *   3) p, the string representation of the permutation.
 *   Returns one if the current decomposition is certified,
 *   otherwise returns zero. */

static PyObject *py2c_factor_update_dobldobl_decomposition
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Updates the decomposition with the given permutation of d elements,
 *   computed in double double precision.
 *   On entry are two integers and one string:
 *   1) d, the number of elements in the permutation;
 *   2) nc, the number of characters in the string;
 *   3) p, the string representation of the permutation.
 *   Returns one if the current decomposition is certified,
 *   otherwise returns zero. */

static PyObject *py2c_factor_update_quaddobl_decomposition
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Updates the decomposition with the given permutation of d elements,
 *   computed in quad double precision.
 *   On entry are two integers and one string:
 *   1) d, the number of elements in the permutation;
 *   2) nc, the number of characters in the string;
 *   3) p, the string representation of the permutation.
 *   Returns one if the current decomposition is certified,
 *   otherwise returns zero. */

static PyObject *py2c_factor_number_of_standard_components
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of irreducible factors in the standard double
 *   precision decomposition of the witness set. */

static PyObject *py2c_factor_number_of_dobldobl_components
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of irreducible factors in the double double
 *   precision decomposition of the witness set. */

static PyObject *py2c_factor_number_of_quaddobl_components
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of irreducible factors in the quad double
 *   precision decomposition of the witness set. */

static PyObject *py2c_factor_witness_points_of_standard_component
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string which represents an irreducible component,
 *   computed in standard double precision.
 *   On entry are two integers:
 *   1) the sum of the degrees of all components;
 *   2) the index of the component. */

static PyObject *py2c_factor_witness_points_of_dobldobl_component
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string which represents an irreducible component,
 *   computed in double double precision.
 *   On entry are two integers:
 *   1) the sum of the degrees of all components;
 *   2) the index of the component. */

static PyObject *py2c_factor_witness_points_of_quaddobl_component
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string which represents an irreducible component,
 *   computed in quad double precision.
 *   On entry are two integers:
 *   1) the sum of the degrees of all components;
 *   2) the index of the component. */

static PyObject *py2c_factor_standard_trace_sum_difference
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the difference between the actual sum at the samples
 *   defined by the labels to the generic points in the factor,
 *   and the trace sum, in standard double precision.
 *   On entry are three integer numbers and one string:
 *   1) d, the number of points in the witness set;
 *   2) k, the dimension of the solution set;
 *   3) nc, the number of characters in the string;
 *   4) ws, the string representing the labels of the witness set. */

static PyObject *py2c_factor_dobldobl_trace_sum_difference
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the difference between the actual sum at the samples
 *   defined by the labels to the generic points in the factor,
 *   and the trace sum, in double double precision.
 *   On entry are three integer numbers and one string:
 *   1) d, the number of points in the witness set;
 *   2) k, the dimension of the solution set;
 *   3) nc, the number of characters in the string;
 *   4) ws, the string representing the labels of the witness set. */

static PyObject *py2c_factor_quaddobl_trace_sum_difference
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the difference between the actual sum at the samples
 *   defined by the labels to the generic points in the factor,
 *   and the trace sum, in quad double precision.
 *   On entry are three integer numbers and one string:
 *   1) d, the number of points in the witness set;
 *   2) k, the dimension of the solution set;
 *   3) nc, the number of characters in the string;
 *   4) ws, the string representing the labels of the witness set. */

static PyObject *py2c_witset_standard_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by an ordinary polynomial system
 *   in standard double precision.
 *   The containers in standard double precision must contain the embedded
 *   polynomial system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 2*nvr doubles for the real and imaginary parts of the
 *   standard double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_dobldobl_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by an ordinary polynomial system
 *   in double double precision.
 *   The containers in double double precision must contain the embedded
 *   polynomial system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 4*nvr doubles for the real and imaginary parts of the
 *   double double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_quaddobl_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by an ordinary polynomial system 
 *   in quad double precision.
 *   The containers in quad double precision must contain the embedded
 *   polynomial system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 8*nvr doubles for the real and imaginary parts of the
 *   quad double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_standard_Laurent_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by a Laurent polynomial system
 *   in standard double precision.
 *   The containers in standard double precision must contain the embedded
 *   Laurent system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 2*nvr doubles for the real and imaginary parts of the
 *   standard double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_dobldobl_Laurent_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by a Laurent polynomial system
 *   in double double precision.
 *   The containers in double double precision must contain the embedded
 *   Laurent system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 4*nvr doubles for the real and imaginary parts of the
 *   double double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_quaddobl_Laurent_membertest
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes the homotopy membership test for a point to belong to
 *   a witness set defined by a Laurent polynomial system
 *   in quad double precision.
 *   The containers in quad double precision must contain the embedded
 *   Laurent system and its corresponding solutions for the witness set
 *   of a positive dimensional solution set.
 *   On entry are the seven parameters, the first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of the point as a list with as
 *   many as 8*nvr doubles for the real and imaginary parts of the
 *   quad double precision coordinates of the test point.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_standard_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by an ordinary polynomial system in standard double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in standard double precision must contain the
 *   embedded system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_dobldobl_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by an ordinary polynomial system in double double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in double double precision must contain the
 *   embedded system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_quaddobl_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by an ordinary polynomial system in quad double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in quad double precision must contain the
 *   embedded system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_standard_Laurent_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by a Laurent polynomial system in standard double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in standard double precision must contain the
 *   embedded Laurent system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_dobldobl_Laurent_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by a Laurent polynomial system in double double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in double double precision must contain the
 *   embedded Laurent system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_witset_quaddobl_Laurent_ismember
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the homotopy membership test for a point to belong to a witness set
 *   defined by a Laurent polynomial system in quad double precision,
 *   where the test point is given as a string in PHCpack format.
 *   The containers in quad double precision must contain the
 *   embedded Laurent system and the corresponding generic points.
 *   On entry are seven parameters.  The first four are integers:
 *   1) vrb, an integer flag (0 or 1) for the verbosity of the test,
 *   2) nvr, the ambient dimension, number of coordinates of the point,
 *   3) dim, the dimension of the witness set,
 *   4) nbc, the number of characters in the string representing the point,
 *   the test point is represented as a solution string in symbolic format,
 *   including the symbols for the variables, before the coordinates;
 *   the next two parameters are two doubles:
 *   5) restol, tolerance on the residual for the valuation of the point,
 *   6) homtol, tolerance on the homotopy membership test for the point;
 *   and the last parameter is a string:
 *   7) tpt, the string representation of a solution which contains
 *   the coordinates of the test point in symbolic format.
 *   On return are three 0/1 integers, to be interpreted as booleans:
 *   1) fail, the failure code of the procedure,
 *   2) onsys, 0 if the evaluation test failed, 1 if success,
 *   3) onset, 0 if not a member of the witness set, 1 if a member. */

static PyObject *py2c_standard_witset_of_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   standard double precision on return contain a witness set for the
 *   hypersurface defined by the ordinary polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of an ordinary polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_witset_of_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   double double precision on return contain a witness set for the
 *   hypersurface defined by the ordinary polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of an ordinary polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_witset_of_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   double double precision on return contain a witness set for the
 *   hypersurface defined by the ordinary polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of an ordinary polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   standard double precision on return contain a witness set for the
 *   hypersurface defined by the Laurent polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of a Laurent polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   double double precision on return contain a witness set for the
 *   hypersurface defined by the Laurent polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of a Laurent polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_witset_of_Laurent_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container in
 *   double double precision on return contain a witness set for the
 *   hypersurface defined by the Laurent polynomial in p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of a Laurent polynomial in several
 *   variables, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_diagonal_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a diagonal homotopy to intersect two solution sets of
 *   dimensions a and b respectively, where a >= b.
 *   The two input parameters are values for a and b.
 *   The systems stored as target and start system in the container,
 *   in standard double precision, define the witness sets for these
 *   two solution sets. */

static PyObject *py2c_dobldobl_diagonal_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a diagonal homotopy to intersect two solution sets of
 *   dimensions a and b respectively, where a >= b.
 *   The two input parameters are values for a and b.
 *   The systems stored as target and start system in the container,
 *   in double double precision, define the witness sets for these
 *   two solution sets. */

static PyObject *py2c_quaddobl_diagonal_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a diagonal homotopy to intersect two solution sets of
 *   dimensions a and b respectively, where a >= b.
 *   The two input parameters are values for a and b.
 *   The systems stored as target and start system in the container,
 *   in quad double precision, define the witness sets for these
 *   two solution sets. */

static PyObject *py2c_standard_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Makes the start solutions to start the cascade homotopy to
 *   intersect two solution sets of dimensions a and b, where a >= b,
 *   in standard double precision.
 *   The dimensions a and b are given as input parameters.
 *   The systems stored as target and start system in the container
 *   define the witness sets for these two solution sets.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_dobldobl_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Makes the start solutions to start the cascade homotopy to
 *   intersect two solution sets of dimensions a and b, where a >= b,
 *   in double double precision.
 *   The dimensions a and b are given as input parameters.
 *   The systems stored as target and start system in the container
 *   define the witness sets for these two solution sets.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_quaddobl_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Makes the start solutions to start the cascade homotopy to
 *   intersect two solution sets of dimensions a and b, where a >= b,
 *   in quad double precision.
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

static PyObject *py2c_diagonal_symbols_doubler 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Doubles the number of symbols in the symbol table to enable the
 *   writing of the target system to string properly when starting the
 *   cascade of a diagonal homotopy in extrinsic coordinates.
 *   On input are three integers, n, d, nc, and one string s.
 *   On input are n, the ambient dimension = #variables before the embedding,
 *   d is the number of slack variables, or the dimension of the first set,
 *   and in s (nc characters) are the symbols for the first witness set.
 *   This function takes the symbols in s and combines those symbols with
 *   those in the current symbol table for the second witness set stored
 *   in the standard systems container.  On return, the symbol table
 *   contains then all symbols to write the top system in the cascade
 *   to start the diagonal homotopy. */

static PyObject *py2c_standard_collapse_diagonal
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Eliminates the extrinsic diagonal for the system and solutions
 *   in the containers for standard doubles.  On input are two integers:
 *   1) k, the current number of slack variables in the embedding;
 *   2) d, the number of slack variables to add to the final embedding.
 *   The system in the container has its diagonal eliminated and is
 *   embedded with k+d slack variables.  The solutions corresponding
 *   to this system are in the solutions container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_collapse_diagonal
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Eliminates the extrinsic diagonal for the system and solutions
 *   in the containers for double doubles.  On input are two integers:
 *   1) k, the current number of slack variables in the embedding;
 *   2) d, the number of slack variables to add to the final embedding.
 *   The system in the container has its diagonal eliminated and is
 *   embedded with k+d slack variables.  The solutions corresponding
 *   to this system are in the solutions container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_collapse_diagonal
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Eliminates the extrinsic diagonal for the system and solutions
 *   in the containers for quad doubles.  On input are two integers:
 *   1) k, the current number of slack variables in the embedding;
 *   2) d, the number of slack variables to add to the final embedding.
 *   The system in the container has its diagonal eliminated and is
 *   embedded with k+d slack variables.  The solutions corresponding
 *   to this system are in the solutions container.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of functions with prototypes in witsols.h starts here. */

static PyObject *py2c_standard_polysys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the standard systems container.  Runs in standard double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_standard_laursys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the standard systems container.  Runs in standard double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_dobldobl_polysys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the dobldobl systems container.  Runs in double double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_dobldobl_laursys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the dobldobl systems container.  Runs in double double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_quaddobl_polysys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the quaddobl systems container.  Runs in quad double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_quaddobl_laursys_solve
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the quaddobl systems container.  Runs in quad double precision.
 *   On input are five integers :
 *   1) nbtasks equals the number of tasks for multitasking,
 *   2) topdim is the top dimension to start the homotopy cascades,
 *   3) filter is a 0 or 1 flag to filter the witness supersets, 
 *   4) factor is a 0 or 1 flag to factor the witness sets,
 *   5) verbose is a flag for intermediate output. */

static PyObject *py2c_copy_standard_polysys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in standard double precision.
 *
 * REQUIRED :
 *   1) py2c_standard_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_copy_standard_laursys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in standard double precision.
 *
 * REQUIRED :
 *   1) py2c_standard_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_copy_dobldobl_polysys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in double double precision.
 *
 * REQUIRED :
 *   1) py2c_dobldobl_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_copy_dobldobl_laursys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in double double precision.
 *
 * REQUIRED :
 *   1) py2c_dobldobl_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_copy_quaddobl_polysys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in quad double precision.
 *
 * REQUIRED :
 *   1) py2c_quaddobl_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_copy_quaddobl_laursys_witset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   There is one integer parameter dim on input,
 *   which represents the dimension of the witness set.
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in quad double precision.
 *
 * REQUIRED :
 *   1) py2c_quaddobl_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

static PyObject *py2c_clear_standard_witsols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in standard double precision. */

static PyObject *py2c_clear_dobldobl_witsols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in double double precision. */

static PyObject *py2c_clear_quaddobl_witsols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in quad double precision. */

/* The wrapping of functions with prototypes in schubert.h starts here. */

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

static PyObject *py2c_schubert_standard_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the Littlewood-Richardson homotopies to resolve a number of
 *   general Schubert intersection conditions on k-planes in n-space,
 *   in standard double precision.
 *   The polynomial system that was solved is in the container for
 *   systems with coefficients in standard double precision and the
 *   corresponding solutions are in the standard solutions container.
 *   On entry are nine integers and two strings, in the following order:
 *   1) n, the ambient dimension, where the k-planes live;
 *   2) k, the dimension of the solution planes;
 *   3) c,the number of intersection conditions;
 *   4) nc, the number of characters in the string brackets;
 *   5) brackets is a string representation of c brackets, where the numbers
 *   in each bracket are separated by spaces;
 *   6) the flag verbose: if 0, then no intermediate output is written,
 *      if 1, then the resolution is dispayed on screen;
 *   7) the flag verify: if 0, then no diagnostic output is written to file,
 *      if 1, then diagnostic output is written to file;
 *   8) the flag minrep: if 0, then all minors are used in the system,
 *      if 1, then a minimal representation of the problem is used;
 *   9) the flag tosquare: if 0, then Gauss-Newton path trackers run,
 *      if 1, then the overdetermined systems are squared;
 *   10) nbchar, the number of characters in the string filename;
 *   11) filename is the name of the output file.
 *   The function returns a tuple of an integer and a string:
 *   0) r is the formal root count as the number of k-planes
 *   for conditions imposed by the brackets for general flags;
 *   1) flags, a string with the coefficients of the general flags. */

static PyObject *py2c_schubert_dobldobl_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the Littlewood-Richardson homotopies to resolve a number of
 *   general Schubert intersection conditions on k-planes in n-space,
 *   in double double precision.
 *   The polynomial system that was solved is in the container for
 *   systems with coefficients in double double precision and the
 *   corresponding solutions are in the dobldobl solutions container.
 *   On entry are nine integers and two strings, in the following order:
 *   1) n, the ambient dimension, where the k-planes live;
 *   2) k, the dimension of the solution planes;
 *   3) c,the number of intersection conditions;
 *   4) nc, the number of characters in the string brackets;
 *   5) brackets is a string representation of c brackets, where the numbers
 *   in each bracket are separated by spaces;
 *   6) the flag verbose: if 0, then no intermediate output is written,
 *      if 1, then the resolution is dispayed on screen;
 *   7) the flag verify: if 0, then no diagnostic output is written to file,
 *      if 1, then diagnostic output is written to file;
 *   8) the flag minrep: if 0, then all minors are used in the system,
 *      if 1, then a minimal representation of the problem is used;
 *   9) the flag tosquare: if 0, then Gauss-Newton path trackers run,
 *      if 1, then the overdetermined systems are squared;
 *   10) nbchar, the number of characters in the string filename;
 *   11) filename is the name of the output file.
 *   The function returns a tuple of an integer and a string:
 *   0) r is the formal root count as the number of k-planes
 *   for conditions imposed by the brackets for general flags;
 *   1) flags, a string with the coefficients of the general flags. */

static PyObject *py2c_schubert_quaddobl_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the Littlewood-Richardson homotopies to resolve a number of
 *   general Schubert intersection conditions on k-planes in n-space,
 *   in quad double precision.
 *   The polynomial system that was solved is in the container for
 *   systems with coefficients in quad double precision and the
 *   corresponding solutions are in the quaddobl solutions container.
 *   On entry are nine integers and two strings, in the following order:
 *   1) n, the ambient dimension, where the k-planes live;
 *   2) k, the dimension of the solution planes;
 *   3) c,the number of intersection conditions;
 *   4) nc, the number of characters in the string brackets;
 *   5) brackets is a string representation of c brackets, where the numbers
 *   in each bracket are separated by spaces;
 *   6) the flag verbose: when 0, then no intermediate output is written,
 *      when 1, then the resolution is dispayed on screen;
 *   7) the flag verify: when 0, then no diagnostic output is written to file,
 *      when 1, then diagnostic output is written to file;
 *   8) the flag minrep: if 0, then all minors are used in the system,
 *      if 1, then a minimal representation of the problem is used;
 *   9) the flag tosquare: if 0, then Gauss-Newton path trackers run,
 *      if 1, then the overdetermined systems are squared;
 *   10) nbchar, the number of characters in the string filename;
 *   11) filename is the name of the output file.
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

/* The wrapping of the functions with prototypes in mapcon.h starts here. */

static PyObject *py2c_mapcon_solve_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Solves the binomial system stored in the Laurent systems container.
 *  There is one input argument, either one or zero.
 *  If one, then only the pure top dimensional solutions are computed.
 *  If zero, then all solution sets are computed.
 *  Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_write_maps
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the maps stored in the container to screen.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_clear_maps
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the maps stored in the container.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_top_dimension
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the top dimension of the maps in the container. */

static PyObject *py2c_mapcon_number_of_maps
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of maps in the container. */

static PyObject *py2c_mapcon_degree_of_map
 ( PyObject *self, PyObject *args );
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

/* The wrapping of functions with prototypes in series.h starts below. */

static PyObject *py2c_standard_Newton_series ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the standard systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Newton_series ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in double double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the double double systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Newton_series ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in quad double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the quad double systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_Newton_power_series
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the standard systems pool the
 *   leading terms of the power series, this function runs Newton's method
 *   to compute power series solutions of the system in the container,
 *   in standard double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the standard systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Newton_power_series
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in double double precision, and in the dobldobl systems pool the
 *   leading terms of the power series, this function runs Newton's method
 *   to compute power series solutions of the system in the container,
 *   in double double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the double double systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Newton_power_series
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in quad double precision, and in the quaddobl systems pool the
 *   leading terms of the power series, this function runs Newton's method
 *   to compute power series solutions of the system in the container,
 *   in quad double precision.  There are four integers on input:
 *   1) the index of the series parameter;
 *   2) the maximal degree of the series;
 *   3) the number of Newton steps to be done on each solution;
 *   4) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The solution series are stored in the quad double systems pool.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_Pade_approximant
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision, followed by the construction of the
 *   Pade approximants, for each solution. There are five integers on input:
 *   1) the index of the series parameter;
 *   2) the degree of the numerator of the Pade approximant;
 *   3) the degree of the denominator of the Pade approximant;
 *   4) the number of Newton steps to be done on each solution;
 *   5) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The Pade approximants are stored in the standard systems pool,
 *   numerators in the odd indexed entries and denominators in the entries
 *   with even index in each system.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Pade_approximant
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in double double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision, followed by the construction of the
 *   Pade approximants, for each solution. There are five integers on input:
 *   1) the index of the series parameter;
 *   2) the degree of the numerator of the Pade approximant;
 *   3) the degree of the denominator of the Pade approximant;
 *   4) the number of Newton steps to be done on each solution;
 *   5) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The Pade approximants are stored in the dobldobl systems pool,
 *   numerators in the odd indexed entries and denominators in the entries
 *   with even index in each system.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Pade_approximant
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in quad double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision, followed by the construction of the
 *   Pade approximants, for each solution. There are five integers on input:
 *   1) the index of the series parameter;
 *   2) the degree of the numerator of the Pade approximant;
 *   3) the degree of the denominator of the Pade approximant;
 *   4) the number of Newton steps to be done on each solution;
 *   5) a 0/1-flag to indicate whether additional diagnostic output needs
 *   to be written to screen.
 *   The Pade approximants are stored in the quaddobl systems pool,
 *   numerators in the odd indexed entries and denominators in the entries
 *   with even index in each system.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of Pade continuation starts here. */

static PyObject *py2c_padcon_set_default_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the default values of the homotopy continuation parameters. */

static PyObject *py2c_padcon_clear_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the allocated space for the parameters. */

static PyObject *py2c_padcon_get_homotopy_continuation_parameter
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the value of the k-th continuation parameter,
 *   if k ranges between 1 and 12.  The integer k is given on entry. */

static PyObject *py2c_padcon_set_homotopy_continuation_gamma
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The gamma constant is the first homotopy continuation parameter.
 *   The gamma is a complex number and it should be given as two
 *   doubles, as its real and imaginary part respectively. */

static PyObject *py2c_padcon_set_homotopy_continuation_parameter
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the value of the k-th continuation parameter to the given value.
 *   The first parameter k is an integer number between 2 and 12.
 *   The second parameter is the value of the k-th parameter,
 *   parsed as a floating point number. */

static PyObject *py2c_padcon_reset_homotopy_continuation_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resets the value of the homotopy continuation parameters
 *   for the step-by-step path trackers.
 *   The first parameter is an integer number, 0, 1, or 2,
 *   respectively for double, double double, or quad double precision. */

static PyObject *py2c_padcon_standard_track
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in standard double precision.
 *   Seven input parameters are expected:
 *   1) the number of characters in the name of the output file;
 *   2) a string which defines the name of the output file,
 *   if the string is empty, then no file is created;
 *   3) a flag to indicate whether the output file is the defined output file
 *   (value 1 of the flag), or whether the file is local (value 0);
 *   4) an integer for the verbose flag, if zero, then no extra
 *   information is written to file or to screen;
 *   5) an integer for the homogenization, if zero, tracking happens in
 *   affine space, if one, then tracking happens in 1-projective space,
 *   if m, for m > 1, then multihomogenization is applied;
 *   6) an integer for the number of variables, 0 if the fifth parameter m
 *   is zero or one;
 *   7) a string with the index representation for the partition of the
 *   set of variables, if the fifth parameter m is larger than one. */

static PyObject *py2c_padcon_dobldobl_track
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in double double precision.
 *   Seven input parameters are expected:
 *   1) the number of characters in the name of the output file;
 *   2) a string which defines the name of the output file,
 *   if the string is empty, then no file is created;
 *   3) a flag to indicate whether the output file is the defined output file
 *   (value 1 of the flag), or whether the file is local (value 0);
 *   4) an integer for the verbose flag, if zero, then no extra
 *   information is written to file or to screen;
 *   5) an integer for the homogenization, if zero, tracking happens in
 *   affine space, if one, then tracking happens in 1-projective space,
 *   if m, for m > 1, then multihomogenization is applied;
 *   6) an integer for the number of variables, 0 if the fifth parameter m
 *   is zero or one;
 *   7) a string with the index representation for the partition of the
 *   set of variables, if the fifth parameter m is larger than one. */

static PyObject *py2c_padcon_quaddobl_track
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in quad double precision.
 *   Seven input parameters are expected:
 *   1) the number of characters in the name of the output file,
 *   2) a string which defines the name of the output file,
 *   if the string is empty, then no file is created;
 *   3) a flag to indicate whether the output file is the defined output file
 *   (value 1 of the flag), or whether the file is local (value 0);
 *   4) an integer for the verbose flag, if zero, then no extra
 *   information is written to file or to screen;
 *   5) an integer for the homogenization, if zero, tracking happens in
 *   affine space, if one, then tracking happens in 1-projective space,
 *   if m, for m > 1, then multihomogenization is applied;
 *   6) an integer for the number of variables, 0 if the fifth parameter m
 *   is zero or one;
 *   7) a string with the index representation for the partition of the
 *   set of variables, if the fifth parameter m is larger than one. */

static PyObject *py2c_padcon_standard_initialize_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in standard double precision,
 *   for the step-by-step Pade continuation.
 *   On entry are two parameters, the verbose flag which is zero or one,
 *   and the homogeneous flag which is zero or one.
 *   If the verbose flag is 1, then extra output will be written.
 *   If the homogeneous flag is 1, tracking happens in projective space. */

static PyObject *py2c_padcon_dobldobl_initialize_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in double double precision,
 *   for the step-by-step Pade continuation.
 *   On entry are two parameters, the verbose flag which is zero or one,
 *   and the homogeneous flag which is zero or one.
 *   If the verbose flag is 1, then extra output will be written.
 *   If the homogeneous flag is 1, tracking happens in projective space. */

static PyObject *py2c_padcon_quaddobl_initialize_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in quad double precision,
 *   for the step-by-step Pade continuation.
 *   On entry is one parameter, the verbose flag which is zero or one,
 *   and the homogeneous flag which is zero or one.
 *   If the verbose flag is 1, then extra output will be written.
 *   If the homogeneous flag is 1, tracking happens in projective space. */

static PyObject *py2c_padcon_standard_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integers: 1) the index for the continuation
 *   parameter in the natural homotopy and 2) the verbose flag.
 *   With the system, defined as target system, and the index
 *   for the continuation parameter, initializes the homotopy in
 *   standard double precision for the step-by-step Pade continuation.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_dobldobl_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integers: 1) the index for the continuation
 *   parameter in the natural homotopy and 2) the verbose flag.
 *   With the system, defined as target system, and the index
 *   for the continuation parameter, initializes the homotopy in
 *   double double precision for the step-by-step Pade continuation.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_quaddobl_initialize_parameter_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integers: 1) the index for the continuation
 *   parameter in the natural homotopy and 2) the verbose flag.
 *   With the system, defined as target system, and the index
 *   for the continuation parameter, initializes the homotopy in
 *   quad double precision for the step-by-step Pade continuation.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_initialize_standard_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Takes the solution with a given index in the solutions container in
 *   standard double precision and initializes the series-Pade tracker.
 *   On entry are two integers: 1) the index of the position of the solution
 *   in the container and 2) the verbose flag, which is zero or one.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_initialize_dobldobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Takes the solution with a given index in the solutions container in
 *   double double precision and initializes the series-Pade tracker.
 *   On entry are two integers: 1) the index of the position of the solution
 *   in the container and 2) the verbose flag, which is zero or one.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_initialize_quaddobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Takes the solution with a given index in the solutions container in
 *   quad double precision and initializes the series-Pade tracker.
 *   On entry are two integers: 1) the index of the position of the solution
 *   in the container and 2) the verbose flag, which is zero or one.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_standard_predict_correct
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in standard double precision.
 *   On entry is one integer, the verbose flag which is zero or one.
 *   On return is the failure code of the predict-correct step:
 *   if zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_dobldobl_predict_correct
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in double double precision.
 *   On entry is one integer, the verbose flag which is zero or one.
 *   On return is the failure code of the predict-correct step:
 *   if zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_quaddobl_predict_correct
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in quad double precision.
 *   On entry is one integer, the verbose flag which is zero or one.
 *   On return is the failure code of the predict-correct step:
 *   if zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_get_standard_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integer parameters: 1) the index of the position of
 *   the solution and 2) the verbose flag, which is zero or one.
 *   Retrieves the current solution and places it at the given position
 *   in the solutions container in standard double precision.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_get_dobldobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integer parameters: 1) the index of the position of
 *   the solution and 2) the verbose flag, which is zero or one.
 *   Retrieves the current solution and places it at the given position
 *   in the solutions container in double double precision.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_get_quaddobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On entry are two integer parameters: 1) the index of the position of
 *   the solution and 2) the verbose flag, which is zero or one.
 *   Retrieves the current solution and places it at the given position
 *   in the solutions container in quad double precision.
 *   If the verbose flag is 1, then extra output will be written. */

static PyObject *py2c_padcon_standard_pole_radius
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the smallest pole radius computed 
 *   by the predictor in standard double precision. */

static PyObject *py2c_padcon_dobldobl_pole_radius
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the smallest pole radius computed
 *   by the predictor in double double precision.
 *   The returned number is the high part of the double double number. */

static PyObject *py2c_padcon_quaddobl_pole_radius
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the smallest pole radius computed
 *   by the predictor in quad double precision.
 *   The returned number is the highest part of the quad double number. */

static PyObject *py2c_padcon_standard_closest_pole
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the complex number representation of the closest pole,
 *   computed by the predictor in standard double precision.
 *   Results are meaningful only if the real part >= 0.0. */

static PyObject *py2c_padcon_dobldobl_closest_pole
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the complex number representation of the closest pole,
 *   computed by the predictor in double double precision.
 *   The returned numbers are the high parts of the double doubles.
 *   Results are meaningful only if the real part >= 0.0.  */

static PyObject *py2c_padcon_quaddobl_closest_pole
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the complex number representation of the closest pole,
 *   computed by the predictor in quad double precision.
 *   The returned numbers are the highest parts of the quad doubles.
 *   Results are meaningful only if the real part >= 0.0. */

static PyObject *py2c_padcon_standard_t_value
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current t value of the path tracker
 *   which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_t_value
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current t value of the path tracker
 *   which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_t_value
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current t value of the path tracker
 *   which runs in quad double precision. */

static PyObject *py2c_padcon_standard_step_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current step size of the path tracker
 *   which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_step_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current step size of the path tracker
 *   which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_step_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current step size of the path tracker
 *   which runs in quad double precision. */

static PyObject *py2c_padcon_standard_series_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current series step size of the path tracker
 *   which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_series_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current series step size of the path tracker
 *   which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_series_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current series step size of the path tracker
 *   which runs in quad double precision. */

static PyObject *py2c_padcon_standard_pole_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current pole step size of the path tracker
 *   which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_pole_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current pole step size of the path tracker
 *   which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_pole_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current pole step size of the path tracker
 *   which runs in quad double precision. */

static PyObject *py2c_padcon_standard_estimated_distance
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the estimated distance to the closest solution by
 *   the path tracker which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_estimated_distance
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the estimated distance to the closest solution by 
 *   the path tracker which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_estimated_distance
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the estimated distance to the closest solution by
 *   the path tracker which runs in quad double precision. */

static PyObject *py2c_padcon_standard_hessian_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current Hessian step size of the path tracker
 *   which runs in standard double precision. */

static PyObject *py2c_padcon_dobldobl_hessian_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current Hessian step size of the path tracker
 *   which runs in double double precision. */

static PyObject *py2c_padcon_quaddobl_hessian_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current Hessian step size of the path tracker
 *   which runs in quad double precision. */

static PyObject *py2c_padcon_standard_series_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx,
 *   of the series computed by the predictor in double precision.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_dobldobl_series_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx, of the
 *   series computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_quaddobl_series_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx, of the
 *   series computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_standard_numerator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the 
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double precision.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_dobldobl_numerator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_quaddobl_numerator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the series
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_standard_denominator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the 
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double precision.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_dobldobl_denominator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_quaddobl_denominator_coefficient
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the series
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles.
 *   The integers leadidx and idx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_standard_pole ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in double precision.
 *   The integers leadidx and poleidx are two input parameters,
 *   the third input integer is the verbose flag. */

static PyObject *py2c_padcon_dobldobl_pole ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in double double precision.
 *   The integers leadidx and poleidx are two input parameters,
 *   the third input integer is the verbose flag.
 *   The returned doubles are the highest parts of the double doubles. */

static PyObject *py2c_padcon_quaddobl_pole ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a tuple: the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in quad double precision.
 *   The integers leadidx and poleidx are two input parameters,
 *   the third input integer is the verbose flag.
 *   The returned doubles are the highest parts of the quad doubles. */

static PyObject *py2c_padcon_clear_standard_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker in double precision. */

static PyObject *py2c_padcon_clear_dobldobl_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker
 *   in double double precision. */

static PyObject *py2c_padcon_clear_quaddobl_data
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker
 *   in quad double precision. */

/* The wrapping of functions with prototypes in syspool.h starts below. */

static PyObject *py2c_syspool_standard_init ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the size of the pool for systems with complex coefficients
 *   in standard double precision, with the value given on input. */

static PyObject *py2c_syspool_dobldobl_init ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the size of the pool for systems with complex coefficients
 *   in double double precision, with the value given on input. */

static PyObject *py2c_syspool_quaddobl_init ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the size of the pool for systems with complex coefficients
 *   in quad double precision, with the value given on input. */

static PyObject *py2c_syspool_standard_size ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the size of the pool for systems in standard double precision. */

static PyObject *py2c_syspool_dobldobl_size ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the size of the pool for systems in double double precision. */

static PyObject *py2c_syspool_quaddobl_size ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the size of the pool for systems in quad double precision. */

static PyObject *py2c_syspool_standard_create
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given the index k (k between one and the size of the pool),
 *   and with a system defined in the standard double system container,
 *   that system is stored as the k-th system in the standard system pool. */

static PyObject *py2c_syspool_dobldobl_create
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given the index k (k between one and the size of the pool),
 *   and with a system defined in the double double system container,
 *   that system is stored as the k-th system in the dobldobl system pool. */

static PyObject *py2c_syspool_quaddobl_create
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given the index k (k between one and the size of the pool),
 *   and with a system defined in the quad double system container,
 *   that system is stored as the k-th system in the quaddobl system pool. */

static PyObject *py2c_syspool_copy_to_standard_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool for systems in standard double
 *   precision to the standard systems container.
 *   The value for k is given as an integer input parameter.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syspool_copy_to_dobldobl_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool for systems in double double
 *   precision to the dobldobl systems container.
 *   The value for k is given as an integer input parameter.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syspool_copy_to_quaddobl_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool for systems in quad double
 *   precision to the quaddobl systems container.
 *   The value for k is given as an integer input parameter.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syspool_standard_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the pool for systems in standard double precision. */

static PyObject *py2c_syspool_dobldobl_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the pool for systems in double double precision. */

static PyObject *py2c_syspool_quaddobl_clear
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the pool for systems in quad double precision. */

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
 *   be copied over from the standard systems container. 
 *   The two other input parameters are two doubles: the real and imaginary
 *   part of the gamma constant.  If the integer parameter equals zero and
 *   if the two input doubles are not both zero, then the input gamma constant
 *   will be used, otherwise, if the two input doubles are zero and the first
 *   integer parameter is zero as well, then a random gamma constant will 
 *   be generated. */

static PyObject *py2c_initialize_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using double double precision arithmetic.
 *   There is one integer number on input to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used.
 *   Before calling this routine the target and start system must
 *   be copied over from the dobldobl systems container.
 *   The two other input parameters are two doubles: the real and imaginary
 *   part of the gamma constant.  If the integer parameter equals zero and
 *   if the two input doubles are not both zero, then the input gamma constant
 *   will be used, otherwise, if the two input doubles are zero and the first
 *   integer parameter is zero as well, then a random gamma constant will 
 *   be generated. */

static PyObject *py2c_initialize_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using quad double precision arithmetic.
 *   There is one integer number on input to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used.
 *   Before calling this routine the target and start system must
 *   be copied over from the quaddobl systems container.
 *   The two other input parameters are two doubles: the real and imaginary
 *   part of the gamma constant.  If the integer parameter equals zero and
 *   if the two input doubles are not both zero, then the input gamma constant
 *   will be used, otherwise, if the two input doubles are zero and the first
 *   integer parameter is zero as well, then a random gamma constant will 
 *   be generated. */

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

static PyObject *py2c_initialize_varbprec_homotopy
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
 * done by algorithmic differentiation in lib2path.h, starts here. */

static PyObject *py2c_ade_newton_d ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The standard systems container contains a valid polynomial system
 *   and the standard solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *            solution container,
 *            if different from zero, then an error happened. */

static PyObject *py2c_ade_newton_dd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation in double
 *   double double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The dobldobl systems container contains a valid polynomial system
 *   and the dobldobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *            solution container,
 *            if different from zero, then an error happened. */

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
 *            solution container,
 *            if different from zero, then an error happened. */

static PyObject *py2c_ade_onepath_d ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the standard solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *            is in the  solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_onepath_dd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation in double
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the dobldobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *            is in the  solution container,
 *            if different from 0, then an error happened. */

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
 *            is in the  solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_manypaths_d ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the standard solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *            are in the solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_manypaths_dd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation in double
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the dobldobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *            are in the solution container,
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
 *            are in the solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_get_default_path_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given the precision (16, 32, or 64), returns the tuple with the
 *   default values for the path parameters.
 *
 * ON ENTRY :
 *   prec     an integer value, 16, 32, or 64 to indicate double,
 *            double double, or quad double precision respectively.
 *
 * ON RETURN :
 *   t        a 14-tuple with the defaults for the path parameters:
 *   t[0]     maximum number of steps along a path;
 *   t[1]     number of points used in the predictor;
 *   t[2]     increase factor of the predictor;
 *   t[3]     decrease factor of the precdictor;
 *   t[4]     maximum step size along a path;
 *   t[5]     maximum step size at the end of a path;
 *   t[6]     minimum step size;
 *   t[7]     tolerance on the residual;
 *   t[8]     tolerance on the corrector update;
 *   t[9]     tolerance on the first correction update;
 *   t[10]    maximum number of iterations of the corrector;
 *   t[11]    tolerance on the corrector;
 *   t[12]    number of steps in the Newton root refiner;
 *   t[13]    tolerance for the final refinement. */

static PyObject *py2c_ade_manypaths_d_pars ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *   All values of the 14 path parameters must be provided,
 *   default values are obtained with py2c_get_default_path_parameters.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the standard solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant;
 *   par00    maximum number of steps along a path;
 *   par01    number of points used in the predictor;
 *   par02    increase factor of the predictor;
 *   par03    decrease factor of the precdictor;
 *   par04    maximum step size along a path;
 *   par05    maximum step size at the end of a path;
 *   par06    minimum step size;
 *   par07    tolerance on the residual;
 *   par08    tolerance on the corrector update;
 *   par09    tolerance on the first correction update;
 *   par10    maximum number of iterations of the corrector;
 *   par11    tolerance on the corrector;
 *   par12    number of steps in the Newton root refiner;
 *   par13    tolerance for the final refinement.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *            are in the solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_manypaths_dd_pars ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation in double
 *   double precision on the data in the systems and solutions container.
 *   All values of the 14 path parameters must be provided,
 *   default values are obtained with py2c_get_default_path_parameters.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the dobldobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant;
 *   par00    maximum number of steps along a path;
 *   par01    number of points used in the predictor;
 *   par02    increase factor of the predictor;
 *   par03    decrease factor of the precdictor;
 *   par04    maximum step size along a path;
 *   par05    maximum step size at the end of a path;
 *   par06    minimum step size;
 *   par07    tolerance on the residual;
 *   par08    tolerance on the corrector update;
 *   par09    tolerance on the first correction update;
 *   par10    maximum number of iterations of the corrector;
 *   par11    tolerance on the corrector;
 *   par12    number of steps in the Newton root refiner;
 *   par13    tolerance for the final refinement.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *            are in the solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_manypaths_qd_pars ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *   All values of the 14 path parameters must be provided,
 *   default values are obtained with py2c_get_default_path_parameters.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant;
 *   par00    maximum number of steps along a path;
 *   par01    number of points used in the predictor;
 *   par02    increase factor of the predictor;
 *   par03    decrease factor of the precdictor;
 *   par04    maximum step size along a path;
 *   par05    maximum step size at the end of a path;
 *   par06    minimum step size;
 *   par07    tolerance on the residual;
 *   par08    tolerance on the corrector update;
 *   par09    tolerance on the first correction update;
 *   par10    maximum number of iterations of the corrector;
 *   par11    tolerance on the corrector;
 *   par12    number of steps in the Newton root refiner;
 *   par13    tolerance for the final refinement.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *            are in the solution container,
 *            if different from 0, then an error happened. */
