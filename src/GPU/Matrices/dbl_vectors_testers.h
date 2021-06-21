// The file dbl_vectors_testers.h specifies test functions
// on vectors of series in double precision.

#ifndef __dbl_vectors_testers_h__
#define __dbl_vectors_testers_h__

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

#endif
