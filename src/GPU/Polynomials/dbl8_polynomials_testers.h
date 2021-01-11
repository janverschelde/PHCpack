/* The dbl8_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in octo double precision. */

#ifndef __dbl8_polynomials_testers_h__
#define __dbl8_polynomials_testers_h__

void dbl8_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo,
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random polynomials and input series.
 *
 * ON ENTRY :
 *   dim         dimension, total number of variables;
 *   nbr         number of terms in the polynomial;
 *   nva         number of variables in each monomial (for products, cyclic);
 *   pwr         highest power of each variable;
 *   deg         truncation degree of the series;
 *   nvr         space for nbr integers;
 *   idx         space for nbr pointers to integers;
 *   exp         space for nbr pointers to integers;
 *   inputhihihi has space for dim arrays of deg+1 doubles;
 *   inputhilohi has space for dim arrays of deg+1 doubles;
 *   inputhihilo has space for dim arrays of deg+1 doubles;
 *   inputhilolo has space for dim arrays of deg+1 doubles;
 *   inputlohihi has space for dim arrays of deg+1 doubles;
 *   inputlolohi has space for dim arrays of deg+1 doubles;
 *   inputlohilo has space for dim arrays of deg+1 doubles;
 *   inputlololo has space for dim arrays of deg+1 doubles;
 *   csthihihi   space for deg+1 doubles;
 *   csthilohi   space for deg+1 doubles;
 *   csthihilo   space for deg+1 doubles;
 *   csthilolo   space for deg+1 doubles;
 *   cstlohihi   space for deg+1 doubles;
 *   cstlolohi   space for deg+1 doubles;
 *   cstlohilo   space for deg+1 doubles;
 *   cstlololo   space for deg+1 doubles;
 *   cffhihihi   space for nbr arrays of deg+1 doubles;
 *   cffhilohi   space for nbr arrays of deg+1 doubles;
 *   cffhihilo   space for nbr arrays of deg+1 doubles;
 *   cffhilolo   space for nbr arrays of deg+1 doubles;
 *   cfflohihi   space for nbr arrays of deg+1 doubles;
 *   cfflolohi   space for nbr arrays of deg+1 doubles;
 *   cfflohilo   space for nbr arrays of deg+1 doubles;
 *   cfflololo   space for nbr arrays of deg+1 doubles;
 *   verbose     if true, then output is written.
 *
 * ON RETURN :
 *   nvr       nvr[k] has the number of variables in monomial k;
 *   idx       idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp       exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputhihihi has the highest doubles of dim input series of degree deg;
 *   inputhilohi has the second highest doubles of dim input series;
 *   inputhihilo has the third highest doubles of dim input series;
 *   inputhilolo has the fourth highest doubles of dim input series;
 *   inputlohihi has the fourth lowest doubles of dim input series;
 *   inputlolohi has the third lowest doubles of dim input series;
 *   inputlohilo has the second lowest doubles of dim input series;
 *   inputlololo has the lowest doubles of dim input series of degree deg;
 *   csthihihi   has the highest doubles of the constant series;
 *   csthilohi   has the second highest doubles of the constant series;
 *   csthihilo   has the third highest doubles of the constant series;
 *   csthilolo   has the fourth highest doubles of the constant series;
 *   cstlohihi   has the fourth lowest doubles of the constant series;
 *   cstlolohi   has the third lowest doubles of the constant series;
 *   cstlohilo   has the second lowest doubles of the constant series;
 *   cstlololo   has the lowest doubles of the constant series;
 *   cffhihihi   cffhihihi[k] has the highest doubles of the coefficient
 *               series of monomial k;
 *   cffhilohi   cffhilohi[k] has the second highest doubles of the
 *               coefficient series of monomial k;
 *   cffhihilo   cffhihilo[k] has the third highest doubles of the coefficient
 *               series of monomial k;
 *   cffhilolo   cffhilolo[k] has the fourth highest doubles of the coefficient
 *               series of monomial k;
 *   cfflohihi   cfflohihi[k] has the fourth lowest doubles of the coefficient
 *               series of monomial k;
 *   cfflolohi   cfflolohi[k] has the third lowest doubles of the coefficient
 *               series of monomial k;
 *   cfflohilo   cfflohilo[k] has the second lowest doubles of the coefficient
 *               series of monomial k;
 *   cfflololo   cfflololo[k] has the lowest doubles of the coefficient series
 *               series of monomial k. */

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
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose,
   bool jobrep=true, int mode=2 );
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
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int main_dbl8_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-120, bool jobrep=true, int mode=2 );
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
 *            this jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int test_dbl8_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode );
/*
 * DESCRIPTION :
 *   For an increasing sequence of degrees,
 *   runs tests in octo double precision.
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
