// The file dbl2_tabs_testers.h specifies test functions on the
// tiled accelerated back substitution in double double precision.

#ifndef __dbl2_tabs_testers_h__
#define __dbl2_tabs_testers_h__

void test_real2_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the computation of its inverse. */

void test_cmplx2_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex upper triangular matrix
 *   to test the computation of its inverse. */

void test_real2_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles and
 *   applies the tiled back substitution to a real random system. */

void test_cmplx2_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles and
 *   applies the tiled back substitution to a complex random system. */

#endif
