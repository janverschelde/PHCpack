// The file dbl_matrices_testers.h specifies test functions on matrices
// and vectors of series in double precision.

#ifndef __dbl_matrices_testers_h__
#define __dbl_matrices_testers_h__

void test_real_inner_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the inner product on random real data. */

void test_cmplx_inner_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the inner product on random complex data. */

void test_real_matrix_vector_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the matrix vector product on random real data. */

void test_cmplx_matrix_vector_product ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a dimension and a degree
 *   and tests the matrix vector product on random complex data. */

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

#endif
