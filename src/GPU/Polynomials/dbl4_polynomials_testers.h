/* The dbl4_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in quad double precision. */

#ifndef __dbl4_polynomials_testers_h__
#define __dbl4_polynomials_testers_h__

void dbl4_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   bool verbose );
/*
 * DESCRIPTION :
 *   Generates random polynomials and input series.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of terms in the polynomial;
 *   nva       number of variables in each monomial (for products, cyclic);
 *   pwr       highest power of each variable;
 *   deg       truncation degree of the series;
 *   nvr       space for nbr integers;
 *   idx       space for nbr pointers to integers;
 *   exp       space for nbr pointers to integers;
 *   inputhihi has space for dim arrays of deg+1 doubles;
 *   inputlohi has space for dim arrays of deg+1 doubles;
 *   inputhilo has space for dim arrays of deg+1 doubles;
 *   inputlolo has space for dim arrays of deg+1 doubles;
 *   csthihi   space for deg+1 doubles;
 *   cstlohi   space for deg+1 doubles;
 *   csthilo   space for deg+1 doubles;
 *   cstlolo   space for deg+1 doubles;
 *   cffhihi   space for nbr arrays of deg+1 doubles;
 *   cfflohi   space for nbr arrays of deg+1 doubles;
 *   cffhilo   space for nbr arrays of deg+1 doubles;
 *   cfflolo   space for nbr arrays of deg+1 doubles;
 *   verbose   if true, then output is written.
 *
 * ON RETURN :
 *   nvr       nvr[k] has the number of variables in monomial k;
 *   idx       idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp       exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputhihi has the highest doubles of dim input series of degree deg;
 *   inputlohi has the second highest doubles of dim input series;
 *   inputhilo has the second lowest doubles of dim input series;
 *   inputlolo has the lowest doubles of dim input series of degree deg;
 *   csthihi   has the highest doubles of the constant series;
 *   cstlohi   has the second highest doubles of the constant series;
 *   csthilo   has the second lowest doubles of the constant series;
 *   cstlolo   has the lowest doubles of the constant series;
 *   cffhihi   cffhihi[k] has the highest doubles of the coefficient series
 *             of monomial k;
 *   cfflohi   cfflohi[k] has the second highest doubles of the coefficient
 *             series of monomial k;
 *   cffhilo   cffhilo[k] has the second lowest doubles of the coefficient
 *             series of monomial k;
 *   cfflolo   cfflolo[k] has the lowest doubles of the coefficient series
 *             of monomial k. */

double dbl4_error_sum
 ( int dim, int deg,
   double **results1hihi_h, double **results1lohi_h, 
   double **results1hilo_h, double **results1lolo_h,
   double **results2hihi_h, double **results2lohi_h,
   double **results2hilo_h, double **results2lolo_h,
   double **resultshihi_d, double **resultslohi_d,
   double **resultshilo_d, double **resultslolo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1hihi_h are the highest doubles computed on the host without jobs;
 *   results1lohi_h are the second highest doubles computed on the host
 *            without jobs;
 *   results1hilo_h are the second lowest doubles computed on the host
 *            without jobs;
 *   results1lolo_h are the lowest doubles computed on the host without jobs;
 *   results2hihi_h are the highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lohi_h are the second highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hilo_h are the second lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lolo_h are the lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   resultshihi_d are the highest doubles computed on the device;
 *   resultslohi_d are the second highest doubles computed on the device;
 *   resultshilo_d are the second lowest doubles computed on the device;
 *   resultslolo_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl4_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose,
   bool jobrep=true );
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
 *            otherwise, the higher the value, the more output;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            this jobs reports needs to be written only once. */

int main_dbl4_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-56, bool jobrep=true );
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
 *            fail if the sum of all errors is larger than tol;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            this jobs reports needs to be written only once. */

#endif
