/* The file next_track.h contains prototypes for the operations to
 * track a solution path with a generator, i.e.: a get_next() method. */

int initialize_standard_homotopy ( int fixed_gamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in standard double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in standard double complex arithmetic. 
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated. */

int initialize_dobldobl_homotopy ( int fixed_gamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in double double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in double double complex arithmetic.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated. */

int initialize_quaddobl_homotopy ( int fixed_gamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in quad double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in quad double complex arithmetic.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated. */

int initialize_multprec_homotopy ( int fixed_gamma, int decimals );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in multiprecision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in multiprecision complex arithmetic with in the working precision
 *   as many decimal places as the value of decimals.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated. */

int initialize_standard_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the standard solution container
 *   and initializes the standard double path tracker with generator. */

int initialize_dobldobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the double double solution container
 *   and initializes the double double path tracker with generator. */

int initialize_quaddobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the quad double solution container
 *   and initializes the double double path tracker with generator. */

int initialize_multprec_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the multiprecision solution container
 *   and initializes the multiprecision path tracker with generator. */

int next_standard_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in standard double precision and replaces the k-th solution in the 
 *   standard solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The standard homotopy has been initialized and the standard path
 *   tracker was initialized with the k-th start solution.  */

int next_dobldobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in double double precision and replaces the k-th solution in the 
 *   double double solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The double double homotopy has been initialized and the double double
 *   path tracker was initialized with the k-th start solution.  */

int next_quaddobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in quad double precision and replaces the k-th solution in the 
 *   quad double solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The quad double homotopy has been initialized and the quad double
 *   path tracker was initialized with the k-th start solution.  */

int next_multprec_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in multiple precision and replaces the k-th solution in the 
 *   multiprecision solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The multiprecision homotopy has been initialized and the multiprecision
 *   path tracker was initialized with the k-th start solution.  */

int clear_standard_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in standard double precision complex arithmetic. */

int clear_dobldobl_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in double double precision complex arithmetic. */

int clear_quaddobl_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in quad double precision complex arithmetic. */

int clear_multprec_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in multiprecision complex arithmetic. */
