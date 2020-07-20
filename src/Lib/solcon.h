/* This file "solcon.h" contains the prototypes of the operations
 * on the solution container in PHCpack.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SOLCON_H__
#define __SOLCON_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int solcon_read_standard_solutions ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file, reads solutions from file, and puts the
 *   solutions in the container; returns 0 if okay, otherwise fail value. */

int solcon_read_standard_solutions_from_file ( int nc, char *filename );
/*
 * DESCRIPTION :
 *   Reads the solutions from the file with name filename, where nc is the
 *   number of characters in the string filename.  The solutions are stored
 *   in the container for the solutions in standard double precision.
 *   Returns 0 if okay, otherwise returns the fail value. */

int solcon_read_dobldobl_solutions ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file, reads solutions from file,
 *   and puts the solutions in the container for double doubles;
 *   returns 0 if okay, otherwise fail value. */

int solcon_read_dobldobl_solutions_from_file ( int nc, char *filename );
/*
 * DESCRIPTION :
 *   Reads the solutions from the file with name filename, where nc is the
 *   number of characters in the string filename.  The solutions are stored
 *   in the container for the solutions in double double precision.
 *   Returns 0 if okay, otherwise returns the fail value. */

int solcon_read_quaddobl_solutions ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file, reads solutions from file,
 *   and puts the solutions in the container for quad doubles;
 *   returns 0 if okay, otherwise fail value. */

int solcon_read_quaddobl_solutions_from_file ( int nc, char *filename );
/*
 * DESCRIPTION :
 *   Reads the solutions from the file with name filename, where nc is the
 *   number of characters in the string filename.  The solutions are stored
 *   in the container for the solutions in quad double precision.
 *   Returns 0 if okay, otherwise returns the fail value. */

int solcon_read_multprec_solutions ( void );
/* 
 * DESCRIPTION :
 *   Prompts the user for a file, reads solutions from file,
 *   and puts the solutions in the multiprecision container;
 *   returns 0 if okay, otherwise fail value. */

int solcon_write_standard_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions in the container to screen. */

int solcon_write_dobldobl_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions in the double double container to screen. */

int solcon_write_quaddobl_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions in the quad double container to screen. */

int solcon_write_multprec_solutions ( void );
/*
 * DESCRIPTION :
 *   Writes the solutions in the multiprecision container to screen. */

int solcon_number_of_standard_solutions ( int *length );
/*
 * DESCRIPTION :
 *   Returns in length the number of the solutions in the container. */

int solcon_number_of_dobldobl_solutions ( int *length );
/*
 * DESCRIPTION :
 *   Returns in length the number of the solutions
 *   in the double double solutions container. */

int solcon_number_of_quaddobl_solutions ( int *length );
/*
 * DESCRIPTION :
 *   Returns in length the number of the solutions
 *   in the quad double solutions container. */

int solcon_number_of_multprec_solutions ( int *length );
/*
 * DESCRIPTION :
 *   Returns in length the number of the solutions
 *   in the multiprecision solutions container. */

int solcon_dimension_of_standard_solutions ( int *dimension );
/*
 * DESCRIPTION :
 *   Returns in dimension the length of the vectors in the container. */

int solcon_dimension_of_dobldobl_solutions ( int *dimension );
/*
 * DESCRIPTION :
 *   Returns in dimension the length of the vectors in 
 *   the container for double double solutions. */

int solcon_dimension_of_quaddobl_solutions ( int *dimension );
/*
 * DESCRIPTION :
 *   Returns in dimension the length of the vectors in 
 *   the container for quad double solutions. */

int solcon_dimension_of_multprec_solutions ( int *dimension );
/*
 * DESCRIPTION :
 *   Returns in dimension the length of the vectors in 
 *   the container for quad double solutions. */

int solcon_retrieve_standard_solution ( int n, int k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the k-th solution n-vector in the container in m and sol.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *
 * ON RETURN :
 *   m        the multiplicity label of the k-th solution;
 *   sol      2*n+5 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 2*n doubles;
 *            3) the last three doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_retrieve_dobldobl_solution ( int n, int k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the k-th solution n-vector in the container in m and sol,
 *   with double double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *
 * ON RETURN :
 *   m        the multiplicity label of the k-th solution;
 *   sol      4*n+10 doubles, with the following meaning:
 *            1) the complex continuation parameter t is the first 4 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 4*n doubles;
 *            3) the last three double doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_retrieve_quaddobl_solution ( int n, int k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the k-th solution n-vector in the container in m and sol,
 *   with quad double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *
 * ON RETURN :
 *   m        the multiplicity label of the k-th solution;
 *   sol      8*n+20 doubles, with the following meaning:
 *            1) the complex continuation parameter t is the first 8 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 8*n doubles;
 *            3) the last three quad doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_retrieve_next_standard_initialize ( void );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current standard solution in the container
 *   to the first solution in the list. */

int solcon_retrieve_next_dobldobl_initialize ( void );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current dobldobl solution in the container
 *   to the first solution in the list. */

int solcon_retrieve_next_quaddobl_initialize ( void );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current quaddobl solution in the container
 *   to the first solution in the list. */

int solcon_retrieve_next_multprec_initialize ( void );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current multprec solution in the container
 *   to the first solution in the list. */

int solcon_retrieve_next_standard_solution
 ( int n, int *k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the next solution n-vector in the container in m and sol.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector.
 *
 * ON RETURN :
 *   k        the position of the solution in the list,
 *            equals zero if the current pointer is null;
 *   m        the multiplicity label of the k-th solution;
 *   sol      2*n+5 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 2*n doubles;
 *            3) the last three doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_retrieve_next_dobldobl_solution
 ( int n, int *k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the next solution n-vector in the container in m and sol,
 *   with double double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector;
 *
 * ON RETURN :
 *   k        the position of the solution in the list,
 *            equals zero if the current pointer is null;
 *   m        the multiplicity label of the k-th solution;
 *   sol      4*n+10 doubles, with the following meaning:
 *            1) the complex continuation parameter t is the first 4 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 4*n doubles;
 *            3) the last three double doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_retrieve_next_quaddobl_solution
 ( int n, int *k, int *m, double *sol );
/* 
 * DESCRIPTION :
 *   Returns the next solution n-vector in the container in m and sol,
 *   with quad double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY :
 *   n        the dimension of the solution vector.
 *
 * ON RETURN :
 *   k        the position of the solution in the list,
 *            equals zero if the current pointer was null;
 *   m        the multiplicity label of the k-th solution;
 *   sol      8*n+20 doubles, with the following meaning:
 *            1) the complex continuation parameter t is the first 8 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 8*n doubles;
 *            3) the last three quad doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_move_current_standard_to_next ( int *cursor );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in standard double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

int solcon_move_current_dobldobl_to_next ( int *cursor );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in double double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

int solcon_move_current_quaddobl_to_next ( int *cursor );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in quad double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

int solcon_move_current_multprec_to_next ( int *cursor );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in multiprecision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

int solcon_length_current_standard_solution_string ( int *cursor, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the current standard double solution in the container,
 *   at the place indicated by the value of the cursor on return.
 *   If this value equals zero, then there is no current solution. */

int solcon_length_current_dobldobl_solution_string ( int *cursor, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the current double double solution in the container,
 *   at the place indicated by the value of the cursor on return.
 *   If this value equals zero, then there is no current solution. */

int solcon_length_current_quaddobl_solution_string ( int *cursor, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the current quad double solution in the container,
 *   at the place indicated by the value of the cursor on return.
 *   If this value equals zero, then there is no current solution. */

int solcon_length_current_multprec_solution_string ( int *cursor, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the current multiprecision solution in the container,
 *   at the place indicated by the value of the cursor on return.
 *   If this value equals zero, then there is no current solution. */

int solcon_write_current_standard_solution_string ( int *k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the current standard double solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0.
 *   The value of k is the place of the solution in the container,
 *   if k equals zero, then there is no current solution. */

int solcon_write_current_dobldobl_solution_string ( int *k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the current solution in the double double solution container to
 *   the string s of n+1 characters.  The last character is \0.
 *   The value of k is the place of the solution in the container,
 *   if k equals zero, then there is no current solution. */

int solcon_write_current_quaddobl_solution_string ( int *k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the current solution in the quad double solution container to
 *   the string s of n+1 characters.  The last character is \0.
 *   The value of k is the place of the solution in the container,
 *   if k equals zero, then there is no current solution. */

int solcon_write_current_multprec_solution_string ( int *k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the current solution in the multiprecision solution container to
 *   the string s of n+1 characters.  The last character is \0.
 *   The value of k is the place of the solution in the container,
 *   if k equals zero, then there is no current solution. */

int solcon_replace_standard_solution ( int n, int k, int m, double *sol );
/*
 * DESCRIPTION :
 *   Replaces the k-th solution n-vector in the container by m and sol,
 *   given in standard double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *   m        the multiplicity label of the k-th solution;
 *   sol      2*n+5 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 2*n doubles;
 *            3) the last three doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_replace_dobldobl_solution ( int n, int k, int m, double *sol );
/*
 * DESCRIPTION :
 *   Replaces the k-th solution n-vector in the container by m and sol,
 *   given in double double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *   m        the multiplicity label of the k-th solution;
 *   sol      4*n+10 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 4*n doubles;
 *            3) the last three double doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_replace_quaddobl_solution ( int n, int k, int m, double *sol );
/*
 * DESCRIPTION :
 *   Replaces the k-th solution n-vector in the container by m and sol,
 *   given in quad double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   k        the position of the solution in the list.
 *   m        the multiplicity label of the k-th solution;
 *   sol      8*n+20 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 4*n doubles;
 *            3) the last three quad doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_append_standard_solution ( int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   Appends a solution n-vector to the container with m and sol.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   m        the multiplicity label of the solution;
 *   sol      2*n+5 doubles, with the following meaning:
 *            1) the complex continuation parameter t is sol[0]+sol[1]*I;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 2*n doubles;
 *            3) the last three doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_append_dobldobl_solution ( int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   Appends a solution n-vector to the container with m and sol,
 *   with double double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   m        the multiplicity label of the solution;
 *   sol      4*n+10 doubles, with the following meaning:
 *            1) the complex continuation parameter t is 
 *               in the first 4 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 4*n doubles;
 *            3) the last three double doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_append_quaddobl_solution ( int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   Appends a solution n-vector to the container with m and sol,
 *   with quad double precision.
 *   The function returns 0 if okay, otherwise it returns the fail value.
 *
 * ON ENTRY:
 *   n        the dimension of the solution vector;
 *   m        the multiplicity label of the solution;
 *   sol      8*n+20 doubles, with the following meaning:
 *            1) the complex continuation parameter t is 
 *               in the first 8 doubles;
 *            2) the real and imaginary parts of the coefficients
 *               of the solution vector are in the next 8*n doubles;
 *            3) the last three quad doubles are respectively
 *               the norm of the last Newton update on the solution;
 *               the inverse of the estimate for the condition number;
 *               the norm of the residual vector. */

int solcon_clear_standard_solutions ( void );
/*
 * DESCRIPTION :
 *   Clears the content of the solution container. */

int solcon_clear_dobldobl_solutions ( void );
/*
 * DESCRIPTION :
 *   Clears the content of the double double solution container. */

int solcon_clear_quaddobl_solutions ( void );
/*
 * DESCRIPTION :
 *   Clears the content of the quad double solution container. */

int solcon_clear_multprec_solutions ( void );
/*
 * DESCRIPTION :
 *   Clears the content of the multiprecision solution container. */

int solcon_open_solution_input_file ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the name of the input file for the solutions and
 *   opens the input file.  All subsequent reading happens from this input. */

int solcon_create_solution_output_file ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the name of the output file for the solutions and
 *   creates the output file, which will be used for all subsequent writing. */

int solcon_scan_solution_banner ( void );
/*
 * DESCRIPTION :
 *   Scans the solution input file for the banner "SOLUTIONS",
 *   this is useful if the solution list is preceded by a system. */

int solcon_write_solution_banner_to_defined_output_file ( void );
/*
 * DESCRIPTION :
 *   Writes the banner "THE SOLUTIONS" to the defined output file. */

int solcon_read_solution_dimensions ( int *len, int *dim );
/*
 * DESCRIPTION :
 *   Reads the dimensions of the solution list from the input file.
 *
 * ON RETURN :
 *   len       length of the solution list, number of solutions;
 *   dim       dimension of the solutions in the list. */

int solcon_write_solution_dimensions ( int len, int dim );
/*
 * DESCRIPTION :
 *   Writes the dimensions of the solution list to file, 
 *   followed by a banner to prepare for more writing of solutions.
 *
 * ON ENTRY :
 *   len       length of the solution list, number of solutions;
 *   dim       dimension of the solutions in the list. */

int solcon_write_solution_dimensions_to_defined_output_file 
      ( int len, int dim );
/*
 * DESCRIPTION :
 *   Writes the dimensions of the solution list to the defined output file,
 *   followed by a banner to prepare for more writing of solutions.
 *
 * ON ENTRY :
 *   len       length of the solution list, number of solutions;
 *   dim       dimension of the solutions in the list. */

int solcon_compute_total_degree_solution
      ( int n, int i, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Computes the i-th solution for the current total degree system.
 *   Returns 1 if the computation failed, otherwise 0 is returned.
 *
 * ON ENTRY :
 *   n         dimension of the solution vector;
 *   i         current solution number.
 *
 * ON RETURN :
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution. */

int solcon_get_linear_product_solution
      ( int n, int i, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Computes the i-th solution for the current linear product system.
 *   Returns 1 if the computation failed, otherwise 0 is returned.
 *
 * ON ENTRY :
 *   n         dimension of the solution vector;
 *   i         current solution number.
 *
 * ON RETURN :
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution. */

int solcon_next_linear_product_solution
      ( int n, int *i, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Computes the next solution for the current linear product system.
 *   Returns 1 if the computation failed, otherwise 0 is returned.
 *
 * ON ENTRY :
 *   n         dimension of the solution vector;
 *   i         current solution counter.
 *
 * ON RETURN :
 *   i         updated solution counter;
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution. */

int solcon_read_next_solution ( int n, int *m, double *sol );
/*
 * DESCRIPTION :
 *   The next solution is read from file.
 *   Returns 1 if reading failed, otherwise 0 is returned.
 *
 * ON ENTRY :
 *   n         dimension of the solution vector.
 *
 * ON RETURN :
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution. */

int solcon_read_next_witness_point ( int k, int n, int *m, double *sol );
/*
 * DESCRIPTION :
 *   The next witness point of set k is read from file.
 *   Returns 1 if reading failed, otherwise 0 is returned.
 *
 * ON ENTRY :
 *   k         index of the witness set, must be 1 or 2;
 *   n         dimension of the solution vector.
 *
 * ON RETURN :
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution. */

int solcon_extrinsic_product
      ( int a, int b, int n1, int n2, double *s1, double *s2,
        int pn, double *ps );
/*
 * DESCRIPTION :
 *   Returns the product of two witness points as start solution to
 *   intersect two witness sets using the extrinsic diagonal homotopy.
 *
 * ON ENTRY :
 *   a         dimension of the first witness set;
 *   b         dimension of the second witness set;
 *   n1        ambient dimension of the first solution;
 *   n2        ambient dimension of the second solution;
 *   s1        witness point of the first set of dimension 2*n1+5;
 *   s2        witness point of the second set of dimension 2*n2+5;
 *   pn        dimension of the product.
 *
 * ON RETURN :
 *   ps        product of the solutions s1 and s2,
 *             of dimension 2*pn+5. */

int add_one_to_double_loop_counters ( int *i, int *j, int n, int m );
/*
 * DESCRIPTION :
 *   Increments j and eventually also i if j reaches the value of m.
 *   Use this routine to invert the double loop: 
 *         for(i=0; i<n; i++) for(j=0; j<m; j++)
 *   into
 *         i=0; j=0;
 *         for(k=0; k<n*m; k++)
 *            done = add_one_to_double_loop_counter(&i,&j,n,m).
 *
 * ON ENTRY :
 *   i        outer counter, ranging between 0 and n-1;
 *   j        inner counter, ranging between 0 and m-1;
 *   n        upper bound for the outer counter;
 *   m        upper bound for the inner counter.
 *
 * ON RETURN :
 *   1 if i has reached n and the loop must stop,
 *   0 otherwise: the double loop continues. */

int solcon_reset_input_file ( int k, int *d, int *n );
/*
 * DESCRIPTION :
 *   Resets the input file to read witness points for set k.
 *
 * ON ENTRY :
 *   k         index of the witness set, must be 1 or 2.
 *
 * ON RETURN :
 *   d         degree of the witness set, equals #solutions;
 *   n         ambient dimension, length of the solution vectors. */

int get_next_start_product
      ( int *i, int *j, int mymonitor,
        int n1, int n2, int dim1, int dim2,
        int deg1, int deg2, int extcd, 
        double *sol1, double *sol2, double *ps );
/*
 * DESCRIPTION :
 *   Reads the next solution of the two witness sets and forms their
 *   product as the next start solution for the diagonal homotopy.
 *
 * ON ENTRY :
 *   i        current value of the outer index in the double loop,
 *            for i ranging between 0 and deg1-1;
 *   j        current value of the inner index in the double loop,
 *            for j ranging between 0 and deg2-1;
 *   monitor  is flag to indicate whether intermediate output or not,
 *            corresponding whether its value is respectively 1 or 0;
 *   n1       ambient dimension of the first witness set;
 *   n2       ambient dimension of the second witness set;
 *   dim1     dimension of the first witness set;
 *   dim2     dimension of the second witness set;
 *   deg1     degree of the first witness set;
 *   deg2     degree of the second witness set;
 *   extcd    dimension of the extrinsic diagonal homotopy;
 *   sol1     workspace for the first solution, of dimension 2*n1+5;
 *   sol2     workspace for the second solution, of dimension 2*n2+5;;
 *   ps       workspace for the product of the first two solutions,
 *            of dimension 2*cd+5.
 *
 * ON RETURN :
 *   sol1     updated next solution of the first witness set;
 *   sol2     updated next solution of the second witness set;
 *   ps       product of the two current start solutions.
 *   
 *   get_next_start_product returns 1 or 0:
 *   1 if some failure occurred or the iteration must stop,
 *   0 if not yet at the end and no failure occurred. */

int solcon_write_next_solution ( int *k, int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   The next solution is written to file.
 *
 * ON ENTRY :
 *   k         current number of solutions already written to file;
 *   n         dimension of the solution vector;
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution.
 *
 * ON RETURN :
 *   k         updated counter on number of solutions written to file. */

int solcon_write_next_solution_to_defined_output_file
      ( int *k, int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   The next solution is written to the defined output file,
 *   as defined by the operation define_output_file in phcpack.h.
 *
 * ON ENTRY :
 *   k         current number of solutions already written to file;
 *   n         dimension of the solution vector;
 *   m         multiplicity flag of the solution vector;
 *   sol       values of the continuation parameter, coordinates,
 *             and diagnostics of the solution.
 *
 * ON RETURN :
 *   k         updated counter on number of solutions written to file. */

int solcon_close_solution_input_file ( int k );
/*
 * DESCRIPTION :
 *   The solution input file is closed for any more reading.
 *   If k equals zero, then the standard input file is closed,
 *   else the file for witness 1 is closed. */

int solcon_close_solution_output_file ( void );
/*
 * DESCRIPTION :
 *   The solution output file is close for any more writing. */

int solcon_length_standard_solution_string ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string
 *   representation of the k-th solution in the container. */

int solcon_length_dobldobl_solution_string ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string
 *   representation of the k-th double double solution in the container. */

int solcon_length_quaddobl_solution_string ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string
 *   representation of the k-th quad double solution in the container. */

int solcon_length_multprec_solution_string ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string
 *   representation of the k-th multiprecision solution in the container. */

int solcon_write_standard_solution_string ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the k-th solution in the solution container to the
 *   string s of n+1 characters.  The last character is \0. */

int solcon_write_dobldobl_solution_string ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the k-th solution in the double double solution container to the
 *   string s of n+1 characters.  The last character is \0. */

int solcon_write_quaddobl_solution_string ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the k-th solution in the quad double solution container to the
 *   string s of n+1 characters.  The last character is \0. */

int solcon_write_multprec_solution_string ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the k-th solution in the multiprecision solution container to the
 *   string s of n+1 characters.  The last character is \0. */

int solcon_length_solution_intro ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the introduction to the k-th solution in the container. */

int solcon_write_solution_intro ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the introduction to the k-th solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0. */

int solcon_length_solution_vector ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the vector in the k-th solution in the container. */

int solcon_write_solution_vector ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the vector in the k-th solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0. */

int solcon_length_solution_diagnostics ( int k, int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of characters in the string representation
 *   of the diagnostics of the k-th solution in the container. */

int solcon_write_solution_diagnostics ( int k, int n, char *s );
/*
 * DESCRIPTION :
 *   Writes the diagnostics of the k-th solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0. */

int solcon_append_standard_solution_string ( int nv, int nc, char *s );
/*
 * DESCRIPTION :
 *   Appends a solution to the container, using the information in the
 *   solution string s of nc+1 characters.  The last character is \0.
 *   The number of variables in the solution equals nv.  */

int solcon_append_dobldobl_solution_string ( int nv, int nc, char *s );
/*
 * DESCRIPTION :
 *   Appends a double double solution to the container, using the data
 *   in the string s of nc+1 characters.  The last character is \0.
 *   The number of variables in the solution equals nv.  */

int solcon_append_quaddobl_solution_string ( int nv, int nc, char *s );
/*
 * DESCRIPTION :
 *   Appends a quad double solution to the container, using the data
 *   in the string s of nc+1 characters.  The last character is \0.
 *   The number of variables in the solution equals nv.  */

int solcon_append_multprec_solution_string ( int nv, int nc, char *s );
/*
 * DESCRIPTION :
 *   Appends a multiprecision solution to the container, using the data
 *   in the string s of nc+1 characters.  The last character is \0.
 *   The number of variables in the solution equals nv.  */

int solcon_replace_solution_string ( int k, int nv, int nc, char *s );
/*
 * DESCRIPTION :
 *   Replaces solution k in the container, using the information in the
 *   solution string s of nc+1 characters.  The last character is \0.
 *   The number of variables in the solution equals nv.  */

int solcon_standard_drop_coordinate_by_index ( int k );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the standard double precision container 
 *   with the same solutions that have their k-th coordinate dropped. */

int solcon_standard_drop_coordinate_by_name ( int nc, char *s );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the standard double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long. */

int solcon_dobldobl_drop_coordinate_by_index ( int k );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the double double precision container 
 *   with the same solutions that have their k-th coordinate dropped. */

int solcon_dobldobl_drop_coordinate_by_name ( int nc, char *s );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the double double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long. */

int solcon_quaddobl_drop_coordinate_by_index ( int k );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the quad double precision container 
 *   with the same solutions that have their k-th coordinate dropped. */

int solcon_quaddobl_drop_coordinate_by_name ( int nc, char *s );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the quad double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long. */

int solcon_standard_set_continuation_parameter ( void );
/*
 * DESCRIPTION :
 *   Sets the value of the continuation parameter to zero for all solutions
 *   in the standard double precision solutions container. */

int solcon_dobldobl_set_continuation_parameter ( void );
/*
 * DESCRIPTION :
 *   Sets the value of the continuation parameter to zero for all solutions
 *   in the double double precision solutions container. */

int solcon_quaddobl_set_continuation_parameter ( void );
/*
 * DESCRIPTION :
 *   Sets the value of the continuation parameter to zero for all solutions
 *   in the quad double precision solutions container. */

int solcon_standard_one_homogenization ( void );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in standard double precision. */

int solcon_dobldobl_one_homogenization ( void );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in double double precision. */

int solcon_quaddobl_one_homogenization ( void );
/*
 * DESCRIPTION :
 *   Add one extra coordinate one to every solution in the container
 *   for solutions in double double precision. */

int solcon_standard_multi_homogenization ( int m );
/*
 * DESCRIPTION :
 *   Add m extra coordinates with value all equal to one to every solution
 *   in the container for solutions in standard double precision. */

int solcon_dobldobl_multi_homogenization ( int m );
/*
 * DESCRIPTION :
 *   Add m extra coordinates with value all equal to one to every solution
 *   in the container for solutions in double double precision. */

int solcon_quaddobl_multi_homogenization ( int m );
/*
 * DESCRIPTION :
 *   Add m extra coordinates with value all equal to one to every solution
 *   in the container for solutions in quad double precision. */

int solcon_standard_one_affinization ( void );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in standard double precision. */

int solcon_dobldobl_one_affinization ( void );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in double double precision. */

int solcon_quaddobl_one_affinization ( void );
/*
 * DESCRIPTION :
 *   Divides every coordinate by the last coordinate of every solution
 *   in the container for solutions in quad double precision. */

int solcon_standard_multi_affinization ( int n, int m, int *z );
/*
 * DESCRIPTION :
 *   For the solutions in the container in standard double precision,
 *   divides coordinates by their corresponding m-homogeneous coordinate.
 *   The number of variables is n, which equals the number of integers
 *   in the index representation z of the partition of the variables. */

int solcon_dobldobl_multi_affinization ( int n, int m, int *z );
/*
 * DESCRIPTION :
 *   For the solutions in the container in double double precision,
 *   divides coordinates by their corresponding m-homogeneous coordinate.
 *   The number of variables is n, which equals the number of integers
 *   in the index representation z of the partition of the variables. */

int solcon_quaddobl_multi_affinization ( int n, int m, int *z );
/*
 * DESCRIPTION :
 *   For the solutions in the container in quad double precision,
 *   divides coordinates by their corresponding m-homogeneous coordinate.
 *   The number of variables is n, which equals the number of integers
 *   in the index representation z of the partition of the variables. */

#endif
