/* The dbl_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in double precision. */

#ifndef __dbl_polynomials_testers_h__
#define __dbl_polynomials_testers_h__

int dbl_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **input, double *cst, double **cff, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random real polynomials and real input series.
 *   Returns 1 if there are duplicates in the support.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables in each monomial (for products, cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   nvr      space for nbr integers;
 *   idx      space for nbr pointers to integers;
 *   exp      space for nbr pointers to integers;
 *   input    space for dim arrays of deg+1 doubles;
 *   cst      space for deg+1 doubles;
 *   cff      space for nbr arrays of deg+1 doubles;
 *   verbose  if true, then output is written.
 *
 * ON RETURN :
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp      exp[k] holds nvr[k] exponents of variables in monomial k;
 *   input    has dim input series of degree deg;
 *   cst      has the coefficients of the constant series;
 *   cff      cff[k] has the coefficient series of monomial k. */

int cmplx_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputre, double **inputim, double *cstre, double *cstim,
   double **cffre, double **cffim, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random complex polynomials and complex input series.
 *   Returns 1 if there are duplicates in the support.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables in each monomial (for products, cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   nvr      space for nbr integers;
 *   idx      space for nbr pointers to integers;
 *   exp      space for nbr pointers to integers;
 *   inputre  space for dim arrays of deg+1 doubles;
 *   inputim  space for dim arrays of deg+1 doubles;
 *   cstre    space for deg+1 doubles;
 *   cstim    space for deg+1 doubles;
 *   cffre    space for nbr arrays of deg+1 doubles;
 *   cffim    space for nbr arrays of deg+1 doubles;
 *   verbose  if true, then output is written.
 *
 * ON RETURN :
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp      exp[k] holds nvr[k] exponents of variables in monomial k;
 *   input    has dim input series of degree deg;
 *   cstre    has the real parts of the coefficients of the constant;
 *   cstim    has the imaginary parts of the coefficients of the constant;
 *   cffre    cffre[k] has the real parts of the coefficient series
 *            of monomial k;
 *   cffim    cffim[k] has the imaginary parts of the coefficient series
 *            of monomial k. */

double dbl_error_sum1
 ( int dim, int deg, double **results_h, double **results_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, for real data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results_h are the results on the host computed without jobs;
 *            computed on the host;
 *   results_d are the results computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double dbl_error_sum
 ( int dim, int deg, double **results1_h, double **results2_h,
   double **results_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, for real data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1_h are the results on the host computed without jobs;
 *   results2_h are the results with convolution and addition jobs,
 *            computed on the host;
 *   results_d are the results computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double cmplx_error_sum
 ( int dim, int deg,
   double **results1re_h, double **results1im_h,
   double **results2re_h, double **results2im_h,
   double **resultsre_d, double **resultsim_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, for complex data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1re_h are the real parts computed on the host without jobs;
 *   results1im_h are the imaginary parts computed on the host without jobs;
 *   results2re_h are the real parts computed on the host,
 *            with convolution and addition jobs;
 *   results2im_h are the imaginary parts computed on the host,
 *            with convolution and addition jobs;
 *   resultsre_d are the real parts computed on the device;
 *   resultsim_d are the imaginary parts computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double cmplx_error_sum1
 ( int dim, int deg,
   double **resultsre_h, double **resultsim_h,
   double **resultsre_d, double **resultsim_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, for complex data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   resultsre_h are the real parts computed on the host without jobs;
 *   resultsim_h are the imaginary parts computed on the host without jobs;
 *   resultsre_d are the real parts computed on the device;
 *   resultsim_d are the imaginary parts computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose,
   bool jobrep=true, int mode=2 );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors if the mode equals 2.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables per monomial (for products and cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written,
 *            otherwise, the higher the value, the more input;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

double test_cmplx_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose,
   bool jobrep=true, int mode=2 );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 *   Returns the sum of all errors if the mode equals 2.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables per monomial (for products and cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written,
 *            otherwise, the higher the value, the more input;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-8, bool jobrep=true, int mode=2 );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests, only if mode is 2.
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
 *   tol      tolerance to decide pass or fail, only for mode equal to 2,
 *            fail if the sum of all errors is larger than tol;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int test_dbl_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode );
/*
 * DESCRIPTION :
 *   For an increasing sequence of degrees,
 *   runs tests in double precision.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables per monomial (for products and cyclic);
 *   pwr      highest power of each variable;
 *   vrblvl   if zero, then no output is written,
 *            otherwise, the higher the value, the more input;
 *   jobrep   if verbose is nonzero and jobrep is true,
 *            then the jobs report is written,
 *            otherwise no jobs report is written.
 *            When running the same problems in many precisions,
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

#endif
