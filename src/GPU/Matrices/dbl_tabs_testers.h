// The file dbl_tabs_testers.h specifies test functions on
// tiled accelerated back substitution in double precision.

#ifndef __dbl_tabs_testers_h__
#define __dbl_tabs_testers_h__

void test_real_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the computation of its inverse. */

void test_cmplx_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random complex upper triangular matrix
 *   to test the computation of its inverse. */

void test_real_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles and
 *   applies the tiled back substitution to a random real system. */

void test_cmplx_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles and
 *   applies the tiled back substitution to a complex random system. */

#endif
