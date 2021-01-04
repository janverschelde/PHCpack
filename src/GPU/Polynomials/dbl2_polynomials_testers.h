/* The dbl2_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in double double precision. */

#ifndef __dbl2_polynomials_testers_h__
#define __dbl2_polynomials_testers_h__

double test_dbl2_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables per monomial (for products and cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written. */

int main_dbl2_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables in each monomial (for products, cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   vrblvl   is the verbose level, if 0 then no output. */

#endif
