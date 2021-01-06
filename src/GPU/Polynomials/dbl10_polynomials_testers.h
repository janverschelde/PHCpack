/* The dbl10_polynomials_testers.h contains the prototypes of functions to
 * test polynomial evaluation and differentiation in deca double precision. */

#ifndef __dbl10_polynomials_testers_h__
#define __dbl10_polynomials_testers_h__

double dbl10_error_sum
 ( int dim, int deg,
   double **results1rtb_h, double **results1rix_h, double **results1rmi_h, 
   double **results1rrg_h, double **results1rpk_h,
   double **results1ltb_h, double **results1lix_h, double **results1lmi_h, 
   double **results1lrg_h, double **results1lpk_h,
   double **results2rtb_h, double **results2rix_h, double **results2rmi_h, 
   double **results2rrg_h, double **results2rpk_h,
   double **results2ltb_h, double **results2lix_h, double **results2lmi_h, 
   double **results2lrg_h, double **results2lpk_h,
   double **resultsrtb_d, double **resultsrix_d, double **resultsrmi_d, 
   double **resultsrrg_d, double **resultsrpk_d,
   double **resultsltb_d, double **resultslix_d, double **resultslmi_d, 
   double **resultslrg_d, double **resultslpk_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1rtb_h are the highest doubles computed on the host without jobs;
 *   results1rix_h are the second highest doubles computed on the host
 *            without jobs;
 *   results1rmi_h are the third highest doubles computed on the host
 *            without jobs;
 *   results1rrg_h are the fourth highest doubles computed on the host
 *            without jobs;
 *   results1rpk_h are the fifth highest doubles computed on the host
 *            without jobs;
 *   results1ltb_h are the fifth lowest doubles computed on the host
 *            without jobs;
 *   results1lix_h are the fourth lowest doubles computed on the host
 *            without jobs;
 *   results1lmi_h are the third lowest doubles computed on the host
 *            without jobs;
 *   results1lrg_h are the second lowest doubles computed on the host
 *            without jobs;
 *   results1lpk_h are the lowest doubles computed on the host without jobs;
 *   results2rtb_h are the highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2rix_h are the second highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2rmi_h are the third highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2rrg_h are the fourth highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2rpk_h are the fifth highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2ltb_h are the fifth lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lix_h are the fourth lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lmi_h are the third lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lrg_h are the second lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lpk_h are the lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   resultsrtb_d are the highest doubles computed on the device;
 *   resultsrix_d are the second highest doubles computed on the device;
 *   resultsrmi_d are the third highest doubles computed on the device;
 *   resultsrrg_d are the fourth highest doubles computed on the device;
 *   resultsrpk_d are the fifth highest doubles computed on the device;
 *   resultsltb_d are the fifth lowest doubles computed on the device;
 *   resultslix_d are the fourth lowest doubles computed on the device;
 *   resultslmi_d are the third lowest doubles computed on the device;
 *   resultslrg_d are the second lowest doubles computed on the device;
 *   resultslpk_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl10_real_polynomial
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

int main_dbl10_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl );
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
 *            if 3 (or higher), then all values are written. */

#endif
