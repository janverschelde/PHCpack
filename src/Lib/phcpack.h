/* file phcpack.h contains prototypes to the operations offered by PHCpack */
/* most BASIC operations in PHCpack : */

int version_string ( int *n, char *s );
/*
 * DESCRIPTION :
 *   Returns in s the string with the current version of PHCpack,
 *   and in n the number of characters filled in s.
 *   The version string takes the form "PHCv2.3.80 released 2013-08-09"
 *   occupying 30 characters.
 *
 * REQUIRED :
 *   Enough space should have been reserved in s for at least n characters. */

int set_seed ( int seed );
/*
 * DESCRIPTION :
 *   Takes the value in seed to initial the seed for the random
 *   number generator. */

int solve_system ( int *root_count );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver of PHCpack on the content of the systems
 *   container.  The solutions on return are in the solution container.
 *   The integer on return is the root count used in the homotopy. */

int solve_Laurent_system ( int *root_count, int silent );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the content of the Laurent systems
 *   container.  The solutions on return are in the solution container.
 *   If silent == 1, then the solver will not write the computed root
 *   counts to screen, otherwise, if silent == 0, the user will see
 *   the computed root counts to screen.
 *   The integer on return is the root count used in the homotopy. */

int mixed_volume ( int *mv );
/*
 * DESCRIPTION :
 *   Computes the mixed volume for the system currently in the systems
 *   container.  The integer in mv on return equals the mixed volume.
 *   The regular mixed-cell configuration is in the cells container. */

int stable_mixed_volume ( int *mv, int *smv );
/*
 * DESCRIPTION :
 *   Computes the mixed volume mv and the stable mixed volume for the 
 *   system currently in the systems container.  The integer in mv on 
 *   return equals the mixed volume.
 *   The regular mixed-cell configuration is in the cells container. */

int deflate ( void );
/*
 * DESCRIPTION :
 *   Applies deflation on the system and solutions in the containers with
 *   standard double precision and default settings of the parameters. */

int Newton_step ( void );
/*
 * DESCRPTION :
 *   Replaces the solutions in the solution container with the results
 *   of one Newton step on the system in the container with the solutions 
 *   in the container on input. */

int dobldobl_Newton_step ( void );
/*
 * DESCRPTION :
 *   Replaces the solutions in the solution container with the results
 *   of one Newton step on the system in the container with the solutions 
 *   in the container on input, using double double arithmetic. */

int quaddobl_Newton_step ( void );
/*
 * DESCRPTION :
 *   Replaces the solutions in the solution container with the results
 *   of one Newton step on the system in the container with the solutions 
 *   in the container on input, using quad double arithmetic. */

int multprec_Newton_step ( int deci );
/*
 * DESCRPTION :
 *   Replaces the solutions in the solution container with the results
 *   of one Newton step on the system in the container with the solutions 
 *   in the container on input, using multiprecision arithmetic.
 *   The input parameter gives the number of decimal places in the
 *   working precision. */

int read_standard_target_system ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   If available on file, also its solutions will be read and stored. */

int read_standard_target_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with standard double precision coefficients
 *   from the file and stores that system into the systems container. */

int read_dobldobl_target_system ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   If available on file, also its solutions will be read and stored.
 *   All data is parsed to double double precision. */

int read_dobldobl_target_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with double double precision coefficients
 *   from the file and stores that system into the systems container. */

int read_quaddobl_target_system ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   If available on file, also its solutions will be read and stored.
 *   All data is parsed to quad double precision. */

int read_quaddobl_target_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with quad double precision coefficients
 *   from the file and stores that system into the systems container. */

int read_multprec_target_system ( int decimals );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   If available on file, also its solutions will be read and stored.
 *   All data is parsed to multiple precision with as many decimal
 *   places as the value of decimals. */

int read_multprec_target_system_from_file
  ( int decimals, int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with multiprecision coefficients
 *   from the file and stores that system into the systems container.
 *   All data is parsed to multiple precision with as many decimal
 *   places as the value of decimals. */

int write_standard_target_system ( void );
/*
 * DESCRIPTION :
 *   Writes the target polynomial system. */

int write_dobldobl_target_system ( void );
/*
 * DESCRIPTION :
 *   Writes the target polynomial system in double double precision. */

int write_quaddobl_target_system ( void );
/*
 * DESCRIPTION :
 *   Writes the target polynomial system in quad double precision. */

int write_multprec_target_system ( void );
/*
 * DESCRIPTION :
 *   Writes the target polynomial system in multiprecision. */

int read_standard_start_system ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file. 
 *   If available, then also start solutions will be read and stored. */

int read_standard_start_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with standard double precision coefficients
 *   and the solutions from the file and stores that system into the 
 *   systems container and the solutions in the solution container. */

int read_dobldobl_start_system ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file. 
 *   If available, then also start solutions will be read and stored.
 *   All data is parsed to double double precision. */

int read_dobldobl_start_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with double double precision coefficients
 *   and the solutions from the file and stores that system into the 
 *   systems container and the solutions in the solution container. */

int read_quaddobl_start_system ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file. 
 *   If available, then also start solutions will be read and stored.
 *   All data is parsed to quad double precision. */

int read_quaddobl_start_system_from_file ( int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with double double precision coefficients
 *   and the solutions from the file and stores that system into the 
 *   systems container and the solutions in the solution container. */

int read_multprec_start_system ( int decimals );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file.
 *   If available on file, also its solutions will be read and stored.
 *   All data is parsed to multiple precision with as many decimal
 *   places as the value of decimals. */

int read_multprec_start_system_from_file
 ( int decimals, int n, char* filename );
/*
 * DESCRIPTION :
 *   Opens the file with name given in the n characters stored in filename,
 *   reads the polynomial system with double double precision coefficients
 *   and the solutions from the file and stores that system into the 
 *   systems container and the solutions in the solution container.
 *   All data is parsed to multiple precision with as many decimal
 *   places as the value of decimals. */

int write_standard_start_system ( void ) ;
/* 
 * DESCRIPTION :
 *   Writes the start polynomial system. */

int write_dobldobl_start_system ( void ) ;
/* 
 * DESCRIPTION :
 *   Writes the start polynomial system in double double precision. */

int write_quaddobl_start_system ( void ) ;
/* 
 * DESCRIPTION :
 *   Writes the start polynomial system in quad double precision. */

int write_multprec_start_system ( void ) ;
/* 
 * DESCRIPTION :
 *   Writes the start polynomial system in multiprecision. */

int write_start_solutions ( void );
/* 
 * DESCRIPTION :
 *   Writes the start solutions. */

int write_dobldobl_start_solutions ( void );
/* 
 * DESCRIPTION :
 *   Writes the start solutions in double double precision. */

int write_quaddobl_start_solutions ( void );
/* 
 * DESCRIPTION :
 *   Writes the start solutions in quad double precision. */

int write_multprec_start_solutions ( void );
/* 
 * DESCRIPTION :
 *   Writes the start solutions in multiprecision. */

int tune_continuation_parameters ( void );
/* 
 * DESCRIPTION :
 *   User can tune the values of the continuation parameters. */

int determine_output_during_continuation ( void );
/* 
 * DESCRIPTION :
 *   User can determine level of output during continuation. */

int retrieve_continuation_parameters ( double *c );
/*
 * DESCRIPTION :
 *   On return are the values of the 34 continuation parameters,
 *   so c must be an array with space for at least 34 doubles. */

int set_continuation_parameters ( double *c );
/*
 * DESCRIPTION :
 *   The values are of the continuation parameters are determined,
 *   based on the 34 values in c. */

int autotune_continuation_parameters
 ( int difficulty_level, int digits_of_precision );
/*
 * DESCRIPTION :
 *   Tunes the continuation parameters based on two parameters:
 *   difficulty_level measures the difficulty level of a path,
 *   with 0 as default, higher values lead to smaller step sizes,
 *   digits_of_precision determines the tolerances along a path. */

int show_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Writes the current values of the continuation parameters to screen. */

int create_homotopy ( void );
/*
 * DESCRIPTION :
 *   Creates a homotopy between the stored target and start system,
 *   for the standard double precision,
 *   using a randomly generated complex number to use as gamma. */

int create_dobldobl_homotopy ( void );
/*
 * DESCRIPTION :
 *   Creates a homotopy between the stored target and start system,
 *   for the double double precision,
 *   using a randomly generated complex number to use as gamma. */

int create_quaddobl_homotopy ( void );
/*
 * DESCRIPTION :
 *   Creates a homotopy between the stored target and start system,
 *   for the quad double precision,
 *   using a randomly generated complex number to use as gamma. */

int create_multprec_homotopy ( void );
/*
 * DESCRIPTION :
 *   Creates a homotopy between the stored target and start system,
 *   for the multiprecision,
 *   using a randomly generated complex number to use as gamma. */

int create_homotopy_with_given_gamma ( double gamma_re, double gamma_im );
/*
 * DESCRIPTION :
 *   Creates a homotopy in standard double precision,
 *   between the stored target and start system,
 *   using the gamma, given by its real and imaginary part. */

int create_dobldobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im );
/*
 * DESCRIPTION :
 *   Creates a homotopy in double double precision,
 *   between the stored target and start system,
 *   using the gamma, given by its real and imaginary part. */

int create_quaddobl_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im );
/*
 * DESCRIPTION :
 *   Creates a homotopy in quad double precision,
 *   between the stored target and start system,
 *   using the gamma, given by its real and imaginary part. */

int create_multprec_homotopy_with_given_gamma
 ( double gamma_re, double gamma_im );
/*
 * DESCRIPTION :
 *   Creates a homotopy in multiprecision,
 *   between the stored target and start system,
 *   using the gamma, given by its real and imaginary part. */

int clear_homotopy ( void );
/*
 * DESCRIPTION :
 *   Clears the homotopy, releasing the allocated memory. */

int clear_dobldobl_homotopy ( void );
/*
 * DESCRIPTION :
 *   Clears the double double homotopy, releasing the allocated memory. */

int clear_quaddobl_homotopy ( void );
/*
 * DESCRIPTION :
 *   Clears the quad double homotopy, releasing the allocated memory. */

int clear_multprec_homotopy ( void );
/*
 * DESCRIPTION :
 *   Clears the multiprecision homotopy, releasing the allocated memory. */

int refine_root ( int n, int *m, double *c );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the system in the container,
 *   to refine the solution given in the input parameters.
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
 *           diagnostics: (err,rco,res) as the last 3 doubles. */

int solve_by_standard_homotopy_continuation ( int number_of_tasks );
/* 
 * DESCRIPTION :
 *   Solves the target system using the start system
 *   and its corresponding start solutions, using as many tasks
 *   as the value of number_of_tasks (if positive).
 *   If number_of_tasks is zero, then no multitasking is used. */

int solve_by_dobldobl_homotopy_continuation ( int number_of_tasks );
/* 
 * DESCRIPTION :
 *   Solves the target system using the start system and its
 *   corresponding start solutions with double double arithmetic,
 *   using as many tasks as the value of number_of_tasks (if positive).
 *   If number_of_tasks is zero, then no multitasking is used. */

int solve_by_quaddobl_homotopy_continuation ( int number_of_tasks );
/* 
 * DESCRIPTION :
 *   Solves the target system using the start system and its
 *   corresponding start solutions with quad double arithmetic,
 *   using as many tasks as the value of number_of_tasks (if positive).
 *   If number_of_tasks is zero, then no multitasking is used. */

int solve_by_multprec_homotopy_continuation ( int decimals );
/* 
 * DESCRIPTION :
 *   Solves the target system using the start system and its
 *   corresponding start solutions with multiprecision arithmetic.
 *   The number of decimal places in the working precision equals
 *   the value given in decimals. */

int write_target_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions of the target system. */

int write_dobldobl_target_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions of the target system in double double precision. */

int write_quaddobl_target_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions of the target system in quad double precision. */

int write_multprec_target_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions of the target system in multiprecision. */

int clear_data ( void );
/* 
 * DESCRIPTION :
 *   Clears the data in PHCpack_Operations. */

int clear_dobldobl_data ( void );
/* 
 * DESCRIPTION :
 *   Clears the double double precision data in PHCpack_Operations. */

int clear_quaddobl_data ( void );
/* 
 * DESCRIPTION :
 *   Clears the quad double precision data in PHCpack_Operations. */

int clear_multprec_data ( void );
/* 
 * DESCRIPTION :
 *   Clears the multiprecision data in PHCpack_Operations. */

int define_output_file ( void );
/*
 * DESCRIPTION :
 *   Asks the user to define an output file. */

int define_output_file_with_string ( int n, char *s );
/*
 * DESCRIPTION :
 *   Opens a file for writing using the string s of n characters.
 *   This file with be the "defined output file". */

int close_output_file ( void );
/*
 * DESCRIPTION :
 *    Closes the defined output file.  */

int write_string_to_defined_output_file ( int n, char *s );
/*
 * DESCRIPTION :
 *   Writes a string s of n characters to the defined output file. */

int write_integers_to_defined_output_file ( int n, int *a );
/*
 * DESCRIPTION :
 *   Writes a sequence a of n integers to the defined output file. */

int write_doubles_to_defined_output_file ( int n, double *a );
/*
 * DESCRIPTION :
 *   Writes a sequence a of n integers to the defined output file. */

/* TRANSFER of data between PHCpack and the containers : */

int copy_target_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies target system to the systems container,
 *   in standard double precision. */

int copy_dobldobl_target_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies target system to the systems container,
 *   in double double precision. */

int copy_quaddobl_target_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies target system to the systems container,
 *   in quad double precision. */

int copy_multprec_target_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies target system to the systems container,
 *   in quad double precision. */

int copy_container_to_target_system ( void );
/* 
 * DESCRIPTION :
 *   Copies system in container to target system,
 *   in standard double precision. */

int copy_dobldobl_container_to_target_system ( void );
/* 
 * DESCRIPTION :
 *   Copies system in container to target system,
 *   in double double precision. */

int copy_quaddobl_container_to_target_system ( void );
/* 
 * DESCRIPTION :
 *   Copies system in container to target system,
 *   in quad double precision. */

int copy_multprec_container_to_target_system ( void );
/* 
 * DESCRIPTION :
 *   Copies system in container to target system,
 *   in quad double precision. */

int copy_start_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start system to the systems container,
 *   in standard double precision. */

int copy_dobldobl_start_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start system to the systems container,
 *   in double double precision. */

int copy_quaddobl_start_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start system to the systems container,
 *   in quad double precision. */

int copy_multprec_start_system_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start system to the systems container,
 *   in quad double precision. */

int copy_container_to_start_system ( void );
/*
 * DESCRIPTION :
 *   Copies system in container to start system,
 *   in standard double precision. */

int copy_dobldobl_container_to_start_system ( void );
/*
 * DESCRIPTION :
 *   Copies system in container to start system,
 *   in double double precision. */

int copy_quaddobl_container_to_start_system ( void );
/*
 * DESCRIPTION :
 *   Copies system in container to start system,
 *   in quad double precision. */

int copy_multprec_container_to_start_system ( void );
/*
 * DESCRIPTION :
 *   Copies system in container to start system,
 *   in quad double precision. */

int copy_target_solutions_to_container ( void );
/* 
 * DESCRIPTION :
 *   Copies target solutions to the solutions container,
 *   in standard double precision. */

int copy_dobldobl_target_solutions_to_container ( void );
/* 
 * DESCRIPTION :
 *   Copies target solutions to the solutions container,
 *   in double double precision. */

int copy_quaddobl_target_solutions_to_container ( void );
/* 
 * DESCRIPTION :
 *   Copies target solutions to the solutions container,
 *   in quad double precision. */

int copy_multprec_target_solutions_to_container ( void );
/* 
 * DESCRIPTION :
 *   Copies target solutions to the solutions container,
 *   in quad double precision. */

int copy_container_to_target_solutions ( void );
/* 
 * DESCRITPION :
 *   Copies solutions in container to target solutions,
 *   in standard double precision. */

int copy_dobldobl_container_to_target_solutions ( void );
/* 
 * DESCRITPION :
 *   Copies solutions in container to target solutions,
 *   in double double precision. */

int copy_quaddobl_container_to_target_solutions ( void );
/* 
 * DESCRITPION :
 *   Copies solutions in container to target solutions,
 *   in quad double precision. */

int copy_multprec_container_to_target_solutions ( void );
/* 
 * DESCRITPION :
 *   Copies solutions in container to target solutions,
 *   in quad double precision. */

int copy_start_solutions_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start solutions to the solutions container,
 *   in standard double precision. */

int copy_dobldobl_start_solutions_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start solutions to the solutions container,
 *   in double double precision. */

int copy_quaddobl_start_solutions_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start solutions to the solutions container,
 *   in quad double precision. */

int copy_multprec_start_solutions_to_container ( void );
/*
 * DESCRIPTION :
 *   Copies start solutions to the solutions container,
 *   in quad double precision. */

int copy_container_to_start_solutions ( void );
/*
 * DESCRIPTION :
 *   Copies solutions in container to start solutions,
 *   in standard double precision. */

int copy_dobldobl_container_to_start_solutions ( void );
/*
 * DESCRIPTION :
 *   Copies solutions in container to start solutions,
 *   in double double precision. */

int copy_quaddobl_container_to_start_solutions ( void );
/*
 * DESCRIPTION :
 *   Copies solutions in container to start solutions,
 *   in quad double precision. */

int copy_multprec_container_to_start_solutions ( void );
/*
 * DESCRIPTION :
 *   Copies solutions in container to start solutions,
 *   in quad double precision. */

/* OPERATIONS on data in the containers : */

int validate_solutions ( void );
/*
 * DESCRIPTION :
 *   Validates solutions in the container,
 *   using the system in the systems container as target. */

int print_system ( void );
/*
 * DESCRIPTION :
 *   Prints the system in the systems container. */
