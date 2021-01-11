/* The dbl10_polynomials_testers.h contains the prototypes of functions to
 * test polynomial evaluation and differentiation in deca double precision. */

#ifndef __dbl10_polynomials_testers_h__
#define __dbl10_polynomials_testers_h__

void dbl10_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double *cstrtb, double *cstrix, double *cstrmi, 
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi, 
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, 
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi, 
   double **cfflrg, double **cfflpk, bool verbose );
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
 *   inputrtb has space for dim arrays of deg+1 doubles;
 *   inputrix has space for dim arrays of deg+1 doubles;
 *   inputrmi has space for dim arrays of deg+1 doubles;
 *   inputrrg has space for dim arrays of deg+1 doubles;
 *   inputrpk has space for dim arrays of deg+1 doubles;
 *   inputltb has space for dim arrays of deg+1 doubles;
 *   inputlix has space for dim arrays of deg+1 doubles;
 *   inputlmi has space for dim arrays of deg+1 doubles;
 *   inputlrg has space for dim arrays of deg+1 doubles;
 *   inputlpk has space for dim arrays of deg+1 doubles;
 *   cstrtb   space for deg+1 doubles;
 *   cstrix   space for deg+1 doubles;
 *   cstrmi   space for deg+1 doubles;
 *   cstrrg   space for deg+1 doubles;
 *   cstrpk   space for deg+1 doubles;
 *   cstltb   space for deg+1 doubles;
 *   cstlix   space for deg+1 doubles;
 *   cstlmi   space for deg+1 doubles;
 *   cstlrg   space for deg+1 doubles;
 *   cstlpk   space for deg+1 doubles;
 *   cffrtb   space for nbr arrays of deg+1 doubles;
 *   cffrix   space for nbr arrays of deg+1 doubles;
 *   cffrmi   space for nbr arrays of deg+1 doubles;
 *   cffrrg   space for nbr arrays of deg+1 doubles;
 *   cffrpk   space for nbr arrays of deg+1 doubles;
 *   cffltb   space for nbr arrays of deg+1 doubles;
 *   cfflix   space for nbr arrays of deg+1 doubles;
 *   cfflmi   space for nbr arrays of deg+1 doubles;
 *   cfflrg   space for nbr arrays of deg+1 doubles;
 *   cfflpk   space for nbr arrays of deg+1 doubles;
 *   verbose  if true, then output is written.
 *
 * ON RETURN :
 *   nvr      nvr[k] has the number of variables in monomial k;
 *   idx      idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp      exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputrtb has the highest doubles of dim input series of degree deg;
 *   inputrix has the second highest doubles of dim input series;
 *   inputrmi has the third highest doubles of dim input series;
 *   inputrrg has the fourth highest doubles of dim input series;
 *   inputrpk has the fifth highest doubles of dim input series;
 *   inputltb has the fifth lowest doubles of dim input series;
 *   inputlix has the fourth lowest doubles of dim input series;
 *   inputlmi has the third lowest doubles of dim input series;
 *   inputlrg has the second lowest doubles of dim input series;
 *   inputlpk has the lowest doubles of dim input series of degree deg;
 *   cstrtb   has the highest doubles of the constant series;
 *   cstrix   has the second highest doubles of the constant series;
 *   cstrmi   has the third highest doubles of the constant series;
 *   cstrrg   has the fourth highest doubles of the constant series;
 *   cstrpk   has the fifth highest doubles of the constant series;
 *   cstltb   has the fifth lowest doubles of the constant series;
 *   cstlix   has the fourth lowest doubles of the constant series;
 *   cstlmi   has the third lowest doubles of the constant series;
 *   cstlrg   has the second lowest doubles of the constant series;
 *   cstlpk   has the lowest doubles of the constant series;
 *   cffrtb   cffrtb[k] has the highest doubles of the coefficient series
 *            of monomial k;
 *   cffrix   cffrix[k] has the second highest doubles of the coefficient
 *            series of monomial k;
 *   cffrmi   cffrmi[k] has the third highest doubles of the coefficient
 *            series of monomial k;
 *   cffrrg   cffrrg[k] has the fourth highest doubles of the coefficient
 *            series of monomial k;
 *   cffrpk   cffrpk[k] has the fifth highest doubles of the coefficient
 *            series of monomial k;
 *   cffltb   cffltb[k] has the fifth lowest doubles of the coefficient
 *            series of monomial k;
 *   cfflix   cfflix[k] has the fourth lowest doubles of the coefficient
 *            series of monomial k;
 *   cfflmi   cfflmi[k] has the third lowest doubles of the coefficient
 *            series of monomial k;
 *   cfflrg   cfflrg[k] has the second lowest doubles of the coefficient
 *            series of monomial k;
 *   cfflpk   cfflpk[k] has the lowest doubles of the coefficient series
 *            of monomial k. */

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

int main_dbl10_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-152, bool jobrep=true, int mode=2 );
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

int test_dbl10_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode );
/*
 * DESCRIPTION :
 *   For an increasing sequence of degrees,
 *   runs tests in deca double precision.
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
