/* The dbl5_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in penta double precision. */

#ifndef __dbl5_polynomials_testers_h__
#define __dbl5_polynomials_testers_h__

double dbl5_error_sum
 ( int dim, int deg,
   double **results1tb_h, double **results1ix_h, double **results1mi_h, 
   double **results1rg_h, double **results1pk_h,
   double **results2tb_h, double **results2ix_h, double **results2mi_h, 
   double **results2rg_h, double **results2pk_h,
   double **resultstb_d, double **resultsix_d, double **resultsmi_d, 
   double **resultsrg_d, double **resultspk_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1tb_h are the highest doubles computed on the host without jobs;
 *   results1ix_h are the second highest doubles computed on the host
 *            without jobs;
 *   results1mi_h are the middle doubles computed on the host without jobs;
 *   results1rg_h are the second lowest doubles computed on the host
 *            without jobs;
 *   results1pk_h are the lowest doubles computed on the host without jobs;
 *   results2tb_h are the highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2ix_h are the second highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2mi_h are the middle doubles computed on the host
 *            with convolution and addition jobs;
 *   results2rg_h are the second lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2pk_h are the lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   resultstb_d are the highest doubles computed on the device;
 *   resultsix_d are the second highest doubles computed on the device;
 *   resultsmi_d are the middle doubles computed on the device;
 *   resultsrg_d are the second lowest doubles computed on the device;
 *   resultspk_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl5_real_polynomial
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

int main_dbl5_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-72 );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in triple double precision.
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
