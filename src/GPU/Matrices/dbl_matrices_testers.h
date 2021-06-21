// The file dbl_matrices_testers.h specifies test functions on matrices
// of series in double precision.

#ifndef __dbl_matrices_testers_h__
#define __dbl_matrices_testers_h__

void test_real_upper_solver ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the backward substitution method. */

void test_cmplx_upper_solver ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex upper triangular matrix
 *   to test the backward substitution method. */

void test_real_lufac ( void );
/*
 * DESCRIPTION :
 *   Generates a random real matrix to test the LU factorization. */

void test_cmplx_lufac ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex matrix to test the LU factorization. */

void test_real_lu_solver ( void );
/*
 * DESCRIPTION :
 *   Generates a random real system to test the LU solver. */

void test_cmplx_lu_solver ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex system to test the LU solver. */

#endif
