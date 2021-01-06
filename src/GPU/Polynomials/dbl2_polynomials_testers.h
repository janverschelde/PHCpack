/* The dbl2_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in double double precision. */

#ifndef __dbl2_polynomials_testers_h__
#define __dbl2_polynomials_testers_h__

double dbl2_error_sum
 ( int dim, int deg,
   double **results1hi_h, double **results1lo_h,
   double **results2hi_h, double **results2lo_h,
   double **resultshi_d, double **resultslo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1hi_h are the high doubles computed on the host without jobs;
 *   results1lo_h are the low doubles computed on the host without jobs;
 *   results2hi_h are the high doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lo_h are the low doubles computed on the host
 *            with convolution and addition jobs;
 *   resultshi_d are the high doubles computed on the device;
 *   resultslo_d are the low doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

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
 *   verbose  if zero, then no output is written,
 *            otherwise, the higher the value, the more output. */

int main_dbl2_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-24 );
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
 *   vrblvl   is the verbose level, if 0 then no output,
 *            otherwise, the higher the value, the more output:
 *            if 1, then the sum of all errors is shown,
 *            if 2, then job counts and timings are listed,
 *            if 3 (or higher), then all values are written;
 *   tol      tolerance to decide pass or fail,
 *            fail if the sum of all errors is larger than tol. */

#endif
