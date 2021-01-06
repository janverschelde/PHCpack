/* The dbl8_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in octo double precision. */

#ifndef __dbl8_polynomials_testers_h__
#define __dbl8_polynomials_testers_h__

double dbl8_error_sum
 ( int dim, int deg,
   double **results1hihihi_h, double **results1hilohi_h, 
   double **results1hihilo_h, double **results1hilolo_h,
   double **results1lohihi_h, double **results1lolohi_h, 
   double **results1lohilo_h, double **results1lololo_h,
   double **results2hihihi_h, double **results2hilohi_h,
   double **results2hihilo_h, double **results2hilolo_h,
   double **results2lohihi_h, double **results2lolohi_h,
   double **results2lohilo_h, double **results2lololo_h,
   double **resultshihihi_d, double **resultshilohi_d,
   double **resultshihilo_d, double **resultshilolo_d,
   double **resultslohihi_d, double **resultslolohi_d,
   double **resultslohilo_d, double **resultslololo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1hihihi_h are the highest doubles computed on the host
 *            without jobs;
 *   results1hilohi_h are the second highest doubles computed on the host
 *            without jobs;
 *   results1hihilo_h are the third highest doubles computed on the host
 *            without jobs;
 *   results1hilolo_h are the fourth highest doubles computed on the host
 *            without jobs;
 *   results1lohihi_h are the fourth lowest doubles computed on the host
 *            without jobs;
 *   results1lolohi_h are the third lowest doubles computed on the host
 *            without jobs;
 *   results1lohilo_h are the second lowest doubles computed on the host
 *            without jobs;
 *   results1lololo_h are the lowest doubles computed on the host
 *            without jobs;
 *   results2hihihi_h are the highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hilohi_h are the second highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hihilo_h are the third highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hilolo_h are the fourth highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lohihi_h are the fourth lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lolohi_h are the third lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lohilo_h are the second lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lololo_h are the lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   resultshihihi_d are the highest doubles computed on the device;
 *   resultshilohi_d are the second highest doubles computed on the device;
 *   resultshihilo_d are the third highest doubles computed on the device;
 *   resultshilolo_d are the fourth highest doubles computed on the device;
 *   resultslohihi_d are the fourth lowest doubles computed on the device;
 *   resultslolohi_d are the third lowest doubles computed on the device;
 *   resultslohilo_d are the second lowest doubles computed on the device;
 *   resultslololo_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl8_real_polynomial
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

int main_dbl8_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-120 );
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
