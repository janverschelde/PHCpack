/* The dbl10_monomials_testers.h contains the prototypes of functions to test
 * monomial evaluation and differentiation in deca double precision. */

#ifndef __dbl10_monomials_testers_h__
#define __dbl10_monomials_testers_h__

double test_dbl10_real ( int dim, int nvr, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written. */

double test_dbl10_complex ( int dim, int nvr, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series.
 *   verbose  if zero, then no output is written. */

int main_dbl10_test
 ( int seed, int dim, int nvr, int pwr, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs two tests on the evaluation and differentiation of a monomial,
 *   on real and complex data.  Returns 0 if all tests are passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   vrblvl   is the verbose level, if 0 then no output. */

#endif
