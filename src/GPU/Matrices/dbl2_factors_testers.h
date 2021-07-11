// The file dbl2_factors_testers.h specifies test functions
// on matrix factorizations in double double precision.

#ifndef __dbl2_factors_testers_h__
#define __dbl2_factors_testers_h__

void test_factors_real2_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization on real data. */

void test_factors_cmplx2_lufac ( void );
/*
 * DESCRIPTION :
 *   Prompts for a dimension and tests the LU factorization
 *   on complex data. */

void test_factors_real2_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on real data. */

void test_factors_cmplx2_houseqr ( void );
/*
 * DESCRIPTION :
 *   Prompts for dimensions and tests the QR decomposition
 *   with Householder matrices on complex data. */

#endif
