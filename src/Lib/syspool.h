/* This file "syspool.h" contains the prototypes of the operations
 * on the systems pool in PHCpack. */

int syspool_standard_initialize ( int n );
/*
 * DESCRIPTION :
 *   Initializes the system pool with n. */

int syspool_standard_size ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the size of the standard systems pool. */

int syspool_dobldobl_size ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the size of the dobldobl systems pool. */

int syspool_quaddobl_size ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the size of the dobldobl systems pool. */

int syspool_standard_read_system ( int k );
/*
 * DESCRIPTION :
 *   Reads the system to create the k-th system in the pool. */

int syspool_standard_write_system ( int k );
/*
 * DESCRIPTION :
 *   Writes the k-th system in the pool. */

int syspool_standard_create ( int k );
/* 
 * DESCRIPTION :
 *   Creates the k-th system in the pool, using the system in the container. */

int syspool_standard_refiner ( int k, int n, int *m, double *c );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the k-th system in the pool,
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

int syspool_copy_to_standard_container ( int k );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool to the standard system container.
 *   Returns the failure code which is zero when all went well. */

int syspool_copy_to_dobldobl_container ( int k );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool to the dobldobl system container.
 *   Returns the failure code which is zero when all went well. */

int syspool_copy_to_quaddobl_container ( int k );
/*
 * DESCRIPTION :
 *   Copies the k-th system in the pool to the quaddobl system container.
 *   Returns the failure code which is zero when all went well. */
