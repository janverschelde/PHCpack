/* The dbl4_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in quad double precision. */

#ifndef __dbl4_polynomials_testers_h__
#define __dbl4_polynomials_testers_h__

int dbl4_make_input
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
 *   Returns 1 if there are duplicates in the support.
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

int cmplx4_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random complex polynomials and complex input series.
 *   Returns 1 if there are duplicates in the support.
 *
 * ON ENTRY :
 *   dim          dimension, total number of variables;
 *   nbr          number of terms in the polynomial;
 *   nva          number of variables in each monomial (for products, cyclic);
 *   pwr          highest power of each variable;
 *   deg          truncation degree of the series;
 *   nvr          space for nbr integers;
 *   idx          space for nbr pointers to integers;
 *   exp          space for nbr pointers to integers;
 *   inputrehihi  space for dim arrays of deg+1 doubles;
 *   inputrelohi  space for dim arrays of deg+1 doubles;
 *   inputrehilo  space for dim arrays of deg+1 doubles;
 *   inputrelolo  space for dim arrays of deg+1 doubles;
 *   inputimhihi  space for dim arrays of deg+1 doubles;
 *   inputimlohi  space for dim arrays of deg+1 doubles;
 *   inputimhilo  space for dim arrays of deg+1 doubles;
 *   inputimlolo  space for dim arrays of deg+1 doubles;
 *   cstrehihi    space for deg+1 doubles;
 *   cstrelohi    space for deg+1 doubles;
 *   cstrehilo    space for deg+1 doubles;
 *   cstrelolo    space for deg+1 doubles;
 *   cstimhihi    space for deg+1 doubles;
 *   cstimlohi    space for deg+1 doubles;
 *   cstimhilo    space for deg+1 doubles;
 *   cstimlolo    space for deg+1 doubles;
 *   cffrehihi    space for nbr arrays of deg+1 doubles;
 *   cffrelohi    space for nbr arrays of deg+1 doubles;
 *   cffrehilo    space for nbr arrays of deg+1 doubles;
 *   cffrelolo    space for nbr arrays of deg+1 doubles;
 *   cffimhihi    space for nbr arrays of deg+1 doubles;
 *   cffimlohi    space for nbr arrays of deg+1 doubles;
 *   cffimhilo    space for nbr arrays of deg+1 doubles;
 *   cffimlolo    space for nbr arrays of deg+1 doubles;
 *   verbose      if true, then output is written.
 *
 * ON RETURN :
 *   nvr          nvr[k] has the number of variables in monomial k;
 *   idx          idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp          exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputrehihi  has the highest doubles of the real parts 
 *                of dim input series of degree deg;
 *   inputrelohi  has the second highest doubles of the real parts 
 *                of dim input series of degree deg;
 *   inputrehilo  has the second lowest doubles of the real parts 
 *                of dim input series of degree deg;
 *   inputrelolo  has the lowest doubles of the real parts
 *                of dim input series of degree deg;
 *   inputimhihi  has the highest doubles of the imaginary parts 
 *                of dim input series of degree deg;
 *   inputimlohi  has the second highest doubles of the imaginary parts 
 *                of dim input series of degree deg;
 *   inputimhilo  has the second lowest doubles of the imaginary parts
 *                of dim input series of degree deg;
 *   inputimlolo  has the lowest doubles of the imaginary parts
 *                of dim input series of degree deg;
 *   cstrehihi    has the highest doubles of the real parts
 *                of the constant series;
 *   cstrelohi    has the second highest doubles of the real parts
 *                of the constant series;
 *   cstrehilo    has the second lowest doubles of the real parts
 *                of the constant series;
 *   cstrelolo    has the lowest doubles of the real parts
 *                of the constant series;
 *   cstimhihi    has the highest doubles of the imaginary parts
 *                of the constant series;
 *   cstimlohi    has the second highest doubles of the imaginary parts
 *                of the constant series;
 *   cstimhilo    has the second lowest doubles of the imaginary parts
 *                of the constant series;
 *   cstimlolo    has the lowest doubles of the imaginary parts
 *                of the constant series;
 *   cffrehihi    cffrehihi[k] has the highest doubles of the real parts
 *                of the coefficient series of monomial k;
 *   cffrelohi    cffrelohi[k] has the second highest doubles of the
 *                real parts of the coefficient series of monomial k;
 *   cffrehilo    cffrehilo[k] has the second lowest doubles of the real parts
 *                of the coefficient series of monomial k;
 *   cffrelolo    cffrelolo[k] has the lowest doubles of the real parts
 *                of the coefficient series of monomial k;
 *   cffimhihi    cffimhihi[k] has the highest doubles of the imaginary parts
 *                of the coefficient series of monomial k;
 *   cffimlohi    cffimlohi[k] has the second highest doubles of the
 *                imaginary parts of the coefficient series of monomial k;
 *   cffimhilo    cffimhilo[k] has the second lowest doubles of the
 *                imaginary parts of the coefficient series of monomial k;
 *   cffimlolo    cffimlo[k] has the lowest doubles of the imaginary parts
 *                of the coefficient series of monomial k. */

double dbl4_error_sum1
 ( int dim, int deg,
   double **resultshihi_h, double **resultslohi_h, 
   double **resultshilo_h, double **resultslolo_h,
   double **resultshihi_d, double **resultslohi_d,
   double **resultshilo_d, double **resultslolo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on real data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   resultshihi_h are the highest doubles computed on the host without jobs;
 *   resultslohi_h are the second highest doubles computed on the host
 *            without jobs;
 *   resultshilo_h are the second lowest doubles computed on the host
 *            without jobs;
 *   resultslolo_h are the lowest doubles computed on the host without jobs;
 *   resultshihi_d are the highest doubles computed on the device;
 *   resultslohi_d are the second highest doubles computed on the device;
 *   resultshilo_d are the second lowest doubles computed on the device;
 *   resultslolo_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

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
 *   with results computed on the device, on real data.
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

double cmplx4_error_sum1
 ( int dim, int deg,
   double **resultsrehihi_h, double **resultsrelohi_h,
   double **resultsrehilo_h, double **resultsrelolo_h,
   double **resultsimhihi_h, double **resultsimlohi_h,
   double **resultsimhilo_h, double **resultsimlolo_h,
   double **resultsrehihi_d, double **resultsrelohi_d,
   double **resultsrehilo_d, double **resultsrelolo_d,
   double **resultsimhihi_d, double **resultsimlohi_d,
   double **resultsimhilo_d, double **resultsimlolo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on complex data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   resultsrehihi_h are the highest doubles of the real parts
 *            computed on the host without jobs;
 *   resultsrelohi_h are the second highest doubles of the real parts
 *            computed on the host without jobs;
 *   resultsrehilo_h are the second lowest doubles of the real parts
 *            computed on the host without jobs;
 *   resultsrelolo_h are the lowest doubles of the real parts
 *            computed on the host without jobs;
 *   resultsimhihi_h are the highest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   resultsimlohi_h are the second highest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   resultsimhilo_h are the second lowest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   resultsimlolo_h are the lowest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   resultsrehihi_d are the highest doubles of the real parts
 *            computed on the device;
 *   resultsrelohi_d are the second highest doubles of the real parts
 *            computed on the device;
 *   resultsrehilo_d are the second lowest doubles of the real parts
 *            computed on the device;
 *   resultsrelolo_d are the lowest doubles of the real parts
 *            computed on the device;
 *   resultsimhihi_d are the highest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimlohi_d are the second highest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimhilo_d are the second lowest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimlolo_d are the lowest doubles of the imaginary parts
 *            computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double cmplx4_error_sum
 ( int dim, int deg,
   double **results1rehihi_h, double **results1relohi_h,
   double **results1rehilo_h, double **results1relolo_h,
   double **results1imhihi_h, double **results1imlohi_h,
   double **results1imhilo_h, double **results1imlolo_h,
   double **results2rehihi_h, double **results2relohi_h,
   double **results2rehilo_h, double **results2relolo_h,
   double **results2imhihi_h, double **results2imlohi_h,
   double **results2imhilo_h, double **results2imlolo_h,
   double **resultsrehihi_d, double **resultsrelohi_d,
   double **resultsrehilo_d, double **resultsrelolo_d,
   double **resultsimhihi_d, double **resultsimlohi_d,
   double **resultsimhilo_d, double **resultsimlolo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on complex data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1rehihi_h are the highest doubles of the real parts
 *            computed on the host without jobs;
 *   results1relohi_h are the second highest doubles of the real parts
 *            computed on the host without jobs;
 *   results1rehilo_h are the second lowest doubles of the real parts
 *            computed on the host without jobs;
 *   results1relolo_h are the lowest doubles of the real parts
 *            computed on the host without jobs;
 *   results1imhihi_h are the highest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   results1imlohi_h are the second highest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   results1imhilo_h are the second lowest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   results1imlolo_h are the lowest doubles of the imaginary parts
 *            computed on the host without jobs;
 *   results2rehihi_h are the highest doubles of the real parts
 *            computed on the host with convolution and addition jobs;
 *   results2relohi_h are the second highest doubles of the real parts
 *            computed on the host with convolution and addition jobs;
 *   results2rehilo_h are the second lowest doubles of the real parts
 *            computed on the host with convolution and addition jobs;
 *   results2relolo_h are the lowest doubles of the real parts
 *            computed on the host with convolution and addition jobs;
 *   results2imhihi_h are the highest doubles of the imaginary parts
 *            computed on the host with convolution and addition jobs;
 *   results2imlohi_h are the second highest doubles of the imaginary parts
 *            computed on the host with convolution and addition jobs;
 *   results2imhilo_h are the second lowest doubles of the imaginary parts
 *            computed on the host with convolution and addition jobs;
 *   results2imlolo_h are the lowest doubles of the imaginary parts
 *            computed on the host with convolution and addition jobs;
 *   resultsrehihi_d are the highest doubles of the real parts
 *            computed on the device;
 *   resultsrelohi_d are the second highest doubles of the real parts
 *            computed on the device;
 *   resultsrehilo_d are the second lowest doubles of the real parts
 *            computed on the device;
 *   resultsrelolo_d are the lowest doubles of the real parts
 *            computed on the device;
 *   resultsimhihi_d are the highest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimlohi_d are the second highest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimhilo_d are the second lowest doubles of the imaginary parts
 *            computed on the device;
 *   resultsimlolo_d are the lowest doubles of the imaginary parts
 *            computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl4_polynomial
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
 *            this jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

double test_cmplx4_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose,
   bool jobrep=true, int mode=2 );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
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

int main_dbl4_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-56, bool jobrep=true, int mode=2 );
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
 *            the jobs reports needs to be written only once;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int test_dbl4_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode );
/*
 * DESCRIPTION :
 *   For an increasing sequence of degrees,
 *   runs tests in quad double precision.
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
