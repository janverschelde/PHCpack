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

void cmplx5_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random complex polynomials and complex input series.
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
 *   inputretb   space for dim arrays of deg+1 doubles;
 *   inputreix   space for dim arrays of deg+1 doubles;
 *   inputremi   space for dim arrays of deg+1 doubles;
 *   inputrerg   space for dim arrays of deg+1 doubles;
 *   inputrepk   space for dim arrays of deg+1 doubles;
 *   inputimtb   space for dim arrays of deg+1 doubles;
 *   inputimix   space for dim arrays of deg+1 doubles;
 *   inputimmi   space for dim arrays of deg+1 doubles;
 *   inputimrg   space for dim arrays of deg+1 doubles;
 *   inputimpk   space for dim arrays of deg+1 doubles;
 *   cstretb     space for deg+1 doubles;
 *   cstreix     space for deg+1 doubles;
 *   cstremi     space for deg+1 doubles;
 *   cstrerg     space for deg+1 doubles;
 *   cstrepk     space for deg+1 doubles;
 *   cstimtb     space for deg+1 doubles;
 *   cstimix     space for deg+1 doubles;
 *   cstimmi     space for deg+1 doubles;
 *   cstimrg     space for deg+1 doubles;
 *   cstimpk     space for deg+1 doubles;
 *   cffretb     space for nbr arrays of deg+1 doubles;
 *   cffreix     space for nbr arrays of deg+1 doubles;
 *   cffremi     space for nbr arrays of deg+1 doubles;
 *   cffrerg     space for nbr arrays of deg+1 doubles;
 *   cffrepk     space for nbr arrays of deg+1 doubles;
 *   cffimtb     space for nbr arrays of deg+1 doubles;
 *   cffimix     space for nbr arrays of deg+1 doubles;
 *   cffimmi     space for nbr arrays of deg+1 doubles;
 *   cffimrg     space for nbr arrays of deg+1 doubles;
 *   cffimpk     space for nbr arrays of deg+1 doubles;
 *   verbose     if true, then output is written.
 *
 * ON RETURN :
 *   nvr         nvr[k] has the number of variables in monomial k;
 *   idx         idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp         exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputretb   has the highest doubles of the real parts 
 *               of dim input series of degree deg;
 *   inputreix   has the second highest doubles of the real parts 
 *               of dim input series of degree deg;
 *   inputremi   has the middle doubles of the real parts 
 *               of dim input series of degree deg;
 *   inputrerg   has the second lowest doubles of the real parts 
 *               of dim input series of degree deg;
 *   inputrepk   has the lowest doubles of the real parts
 *               of dim input series of degree deg;
 *   inputimtb   has the highest doubles of the imaginary parts 
 *               of dim input series of degree deg;
 *   inputimix   has the second highest doubles of the imaginary parts 
 *               of dim input series of degree deg;
 *   inputimmi   has the middle doubles of the imaginary parts 
 *               of dim input series of degree deg;
 *   inputimrg   has the second lowest doubles of the imaginary parts
 *               of dim input series of degree deg;
 *   inputimpk   has the lowest doubles of the imaginary parts
 *               of dim input series of degree deg;
 *   cstretb     has the highest doubles of the real parts
 *               of the constant series;
 *   cstreix     has the second highest doubles of the real parts
 *               of the constant series;
 *   cstremi     has the middle doubles of the real parts
 *               of the constant series;
 *   cstrerg     has the second lowest doubles of the real parts
 *               of the constant series;
 *   cstrepk     has the lowest doubles of the real parts
 *               of the constant series;
 *   cstimtb     has the highest doubles of the imaginary parts
 *               of the constant series;
 *   cstimix     has the second highest doubles of the imaginary parts
 *               of the constant series;
 *   cstimmi     has the middle doubles of the imaginary parts
 *               of the constant series;
 *   cstimrg     has the second lowest doubles of the imaginary parts
 *               of the constant series;
 *   cstimpk     has the lowest doubles of the imaginary parts
 *               of the constant series;
 *   cffretb     cffretb[k] has the highest doubles of the real parts
 *               of the coefficient series of monomial k;
 *   cffreix     cffreix[k] has the second highest doubles of the
 *               real parts of the coefficient series of monomial k;
 *   cffremi     cffremi[k] has the middle doubles of the
 *               real parts of the coefficient series of monomial k;
 *   cffrerg     cffrerg[k] has the second lowest doubles of the real parts
 *               of the coefficient series of monomial k;
 *   cffrepk     cffrepk[k] has the lowest doubles of the real parts
 *               of the coefficient series of monomial k;
 *   cffimtb     cffimtb[k] has the highest doubles of the imaginary parts
 *               of the coefficient series of monomial k;
 *   cffimix     cffimix[k] has the second highest doubles of the
 *               imaginary parts of the coefficient series of monomial k;
 *   cffimmi     cffimmi[k] has the middle doubles of the
 *               imaginary parts of the coefficient series of monomial k;
 *   cffimrg     cffimrg[k] has the second lowest doubles of the
 *               imaginary parts of the coefficient series of monomial k;
 *   cffimpk     cffimlo[k] has the lowest doubles of the imaginary parts
 *               of the coefficient series of monomial k. */

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

double test_cmplx5_real_polynomial
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

int main_dbl5_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol=1.0e-72, bool jobrep=true, int mode=2 );
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

int test_dbl5_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode );
/*
 * DESCRIPTION :
 *   For an increasing sequence of degrees,
 *   runs tests in penta double precision.
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
