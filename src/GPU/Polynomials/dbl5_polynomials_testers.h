/* The dbl5_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in penta double precision. */

#ifndef __dbl5_polynomials_testers_h__
#define __dbl5_polynomials_testers_h__

void dbl5_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double *csttb, double *cstix, double *cstmi, 
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, 
   double **cffrg, double **cffpk, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random polynomials and input series.
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
 *   inputtb  space for dim arrays of deg+1 doubles;
 *   inputix  space for dim arrays of deg+1 doubles;
 *   inputmi  space for dim arrays of deg+1 doubles;
 *   inputrg  space for dim arrays of deg+1 doubles;
 *   inputpk  space for dim arrays of deg+1 doubles;
 *   csttb    space for deg+1 doubles;
 *   cstix    space for deg+1 doubles;
 *   cstmi    space for deg+1 doubles;
 *   cstrg    space for deg+1 doubles;
 *   cstpk    space for deg+1 doubles;
 *   cfftb    space for nbr arrays of deg+1 doubles;
 *   cffix    space for nbr arrays of deg+1 doubles;
 *   cffmi    space for nbr arrays of deg+1 doubles;
 *   cffrg    space for nbr arrays of deg+1 doubles;
 *   cffpk    space for nbr arrays of deg+1 doubles;
 *   verbose  if true, then output is written.
 *
 * ON RETURN :
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp      exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputtb  has the highest doubles of dim input series of degree deg;
 *   inputix  has the second highest doubles of dim input series;
 *   inputmi  has the middle doubles of dim input series of degree deg;
 *   inputrg  has the second lowest doubles of dim input series;
 *   inputpk  has the lowest doubles of dim input series of degree deg;
 *   csttb    has the highest doubles of the constant series;
 *   cstix    has the second highest doubles of the constant series;
 *   cstmi    has the middle doubles of the constant series;
 *   cstrg    has the second lowest doubles of the constant series;
 *   cstpk    has the lowest doubles of the constant series;
 *   cfftb    cfftb[k] has the highest doubles of the coefficient series
 *            of monomial k;
 *   cffix    cffix[k] has the second highest doubles of the coefficient
 *            series of monomial k;
 *   cffmi    cffmi[k] has the middle doubles of the coefficient series
 *            of monomial k;
 *   cffrg    cffrg[k] has the second lowest doubles of the coefficient
 *            series of monomial k;
 *   cffpk    cffpk[k] has the lowest doubles of the coefficient series
 *            of monomial k. */

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
