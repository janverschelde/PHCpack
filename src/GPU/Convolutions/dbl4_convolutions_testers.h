/* The file dbl4_convolutions_testers.h contains the prototypes of
 * functions to test the product of two series in quad double precision. */

#ifndef __dbl4_convolutions_testers_h__
#define __dbl4_convolutions_testers_h__

double test_dbl4_real ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for real coefficients.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

double test_dbl4_real_random ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies two random power series,
 *   truncated to degree deg, for real coefficients.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

double test_dbl4_complex ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for complex coefficients.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

double test_dbl4_complex_random ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies two random power series,
 *   truncated to degree deg, for complex coefficients.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

double test_dbl4_real_exponential ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random x in [-1,+1], for real coefficients
 *   of a series of degree truncated to deg.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

double test_dbl4_complex_exponential ( int deg, int verbose );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random complex number on the unit circle,
 *   for series of degree truncated to deg.
 *   If verbose equals zero, then no output is written.
 *   Returns the sum of all errors. */

int main_dbl4_test ( int seed, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs six tests on convolutions in quad double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed of the random number generators,
 *            if 0, then the current time will be used as seed;
 *   deg      degree of all series;
 *   vrblvl   is the verbose level:
 *            if 0, then there is no output,
 *            if 1, then only the pass/fail conclusions are written,
 *            if 2, then the numerical results are shown. */

#endif
