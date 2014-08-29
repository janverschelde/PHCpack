/* file unisolvers.h contains prototypes to the univariate root finders,
 * wrapped through the Ada code unisolve of PHCpack */

int solve_with_standard_doubles ( int max, double eps, int *nit );
/*
 * DESCRIPTION :
 *   Calls the method of Durand-Kerner as implemented in PHCpack
 *   in standard double precision.
 *
 * REQUIRED :
 *   The standard polynomial systems container must have as system
 *   one polynomial in one variable, before calling this function.
 *   On return, approximations for the roots are in the standard 
 *   solutions container.
 *
 * ON ENTRY :
 *   max      maximum number of iterations;
 *   eps      accuracy requirement.
 *
 * ON RETURN :
 *   nit      number of iterations performed by Durand-Kerner. */

int solve_with_double_doubles ( int max, double eps, int *nit );
/*
 * DESCRIPTION :
 *   Calls the method of Durand-Kerner as implemented in PHCpack
 *   in double double precision.
 *
 * REQUIRED :
 *   The dobldobl polynomial systems container must have as system
 *   one polynomial in one variable, before calling this function.
 *   On return, approximations for the roots are in the dobldobl
 *   solutions container.
 *
 * ON ENTRY :
 *   max      maximum number of iterations;
 *   eps      accuracy requirement.
 *
 * ON RETURN :
 *   nit      number of iterations performed by Durand-Kerner. */

int solve_with_quad_doubles ( int max, double eps, int *nit );
/*
 * DESCRIPTION :
 *   Calls the method of Durand-Kerner as implemented in PHCpack
 *   in quad double precision.
 *
 * REQUIRED :
 *   The quaddobl polynomial systems container must have as system
 *   one polynomial in one variable, before calling this function.
 *   On return, approximations for the roots are in the quaddobl
 *   solutions container.
 *
 * ON ENTRY :
 *   max      maximum number of iterations;
 *   eps      accuracy requirement.
 *
 * ON RETURN :
 *   nit      number of iterations performed by Durand-Kerner. */

int solve_with_multiprecision ( int dcp, int max, double eps, int *nit );
/*
 * DESCRIPTION :
 *   Calls the method of Durand-Kerner as implemented in PHCpack
 *   in quad double precision.
 *
 * REQUIRED :
 *   The multprec polynomial systems container must have as system
 *   one polynomial in one variable, before calling this function.
 *   On return, approximations for the roots are in the multprec
 *   solutions container.
 *
 * ON ENTRY :
 *   dcp      number of decimal places in the working precision;
 *   max      maximum number of iterations;
 *   eps      accuracy requirement.
 *
 * ON RETURN :
 *   nit      number of iterations performed by Durand-Kerner. */
