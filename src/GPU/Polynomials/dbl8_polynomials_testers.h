/* The dbl8_polynomials_testers.h contains the prototypes of functions to test
 * polynomial evaluation and differentiation in octo double precision. */

#ifndef __dbl8_polynomials_testers_h__
#define __dbl8_polynomials_testers_h__

int dbl8_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random polynomials and input series for real data.
 *   Returns 1 if there are duplicates in the support.
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
 *   inputlohihi has space for dim arrays of deg+1 doubles;
 *   inputhilohi has space for dim arrays of deg+1 doubles;
 *   inputlolohi has space for dim arrays of deg+1 doubles;
 *   inputhihilo has space for dim arrays of deg+1 doubles;
 *   inputlohilo has space for dim arrays of deg+1 doubles;
 *   inputhilolo has space for dim arrays of deg+1 doubles;
 *   inputlololo has space for dim arrays of deg+1 doubles;
 *   csthihihi   space for deg+1 doubles;
 *   cstlohihi   space for deg+1 doubles;
 *   csthilohi   space for deg+1 doubles;
 *   cstlolohi   space for deg+1 doubles;
 *   csthihilo   space for deg+1 doubles;
 *   cstlohilo   space for deg+1 doubles;
 *   csthilolo   space for deg+1 doubles;
 *   cstlololo   space for deg+1 doubles;
 *   cffhihihi   space for nbr arrays of deg+1 doubles;
 *   cfflohihi   space for nbr arrays of deg+1 doubles;
 *   cffhilohi   space for nbr arrays of deg+1 doubles;
 *   cfflolohi   space for nbr arrays of deg+1 doubles;
 *   cffhihilo   space for nbr arrays of deg+1 doubles;
 *   cfflohilo   space for nbr arrays of deg+1 doubles;
 *   cffhilolo   space for nbr arrays of deg+1 doubles;
 *   cfflololo   space for nbr arrays of deg+1 doubles;
 *   verbose     if true, then output is written.
 *
 * ON RETURN :
 *   nvr       nvr[k] has the number of variables in monomial k;
 *   idx       idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp       exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputhihihi has the highest doubles of dim input series of degree deg;
 *   inputlohihi has the second highest doubles of dim input series;
 *   inputhilohi has the third highest doubles of dim input series;
 *   inputlolohi has the fourth highest doubles of dim input series;
 *   inputhihilo has the fourth lowest doubles of dim input series;
 *   inputlohilo has the third lowest doubles of dim input series;
 *   inputhilolo has the second lowest doubles of dim input series;
 *   inputlololo has the lowest doubles of dim input series of degree deg;
 *   csthihihi   has the highest doubles of the constant series;
 *   cstlohihi   has the second highest doubles of the constant series;
 *   csthilohi   has the third highest doubles of the constant series;
 *   cstlolohi   has the fourth highest doubles of the constant series;
 *   csthihilo   has the fourth lowest doubles of the constant series;
 *   cstlohilo   has the third lowest doubles of the constant series;
 *   csthilolo   has the second lowest doubles of the constant series;
 *   cstlololo   has the lowest doubles of the constant series;
 *   cffhihihi   cffhihihi[k] has the highest doubles of the coefficient
 *               series of monomial k;
 *   cfflohihi   cfflohihi[k] has the second highest doubles of the
 *               coefficient series of monomial k;
 *   cffhilohi   cffhilohi[k] has the third highest doubles of the
 *               coefficient series of monomial k;
 *   cfflolohi   cfflolohi[k] has the fourth highest doubles of the
 *               coefficient series of monomial k;
 *   cffhihilo   cffhihilo[k] has the fourth lowest doubles of the
 *               coefficient series of monomial k;
 *   cfflohilo   cfflohilo[k] has the third lowest doubles of the
 *               coefficient series of monomial k;
 *   cffhilolo   cffhilolo[k] has the second lowest doubles of the
 *               coefficient series of monomial k;
 *   cfflololo   cfflololo[k] has the lowest doubles of the coefficient
 *               series of monomial k. */

int cmplx8_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo, bool verbose );
/*
 * DESCRIPTION :
 *   Generates random polynomials and input series for complex data.
 *   Returns 1 if there are duplicates in the support.
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
 *   inputrehihihi has space for dim arrays of deg+1 doubles;
 *   inputrelohihi has space for dim arrays of deg+1 doubles;
 *   inputrehilohi has space for dim arrays of deg+1 doubles;
 *   inputrelolohi has space for dim arrays of deg+1 doubles;
 *   inputrehihilo has space for dim arrays of deg+1 doubles;
 *   inputrelohilo has space for dim arrays of deg+1 doubles;
 *   inputrehilolo has space for dim arrays of deg+1 doubles;
 *   inputrelololo has space for dim arrays of deg+1 doubles;
 *   inputimhihihi has space for dim arrays of deg+1 doubles;
 *   inputimlohihi has space for dim arrays of deg+1 doubles;
 *   inputimhilohi has space for dim arrays of deg+1 doubles;
 *   inputimlolohi has space for dim arrays of deg+1 doubles;
 *   inputimhihilo has space for dim arrays of deg+1 doubles;
 *   inputimlohilo has space for dim arrays of deg+1 doubles;
 *   inputimhilolo has space for dim arrays of deg+1 doubles;
 *   inputimlololo has space for dim arrays of deg+1 doubles;
 *   cstrehihihi has space for deg+1 doubles;
 *   cstrelohihi has space for deg+1 doubles;
 *   cstrehilohi has space for deg+1 doubles;
 *   cstrelolohi has space for deg+1 doubles;
 *   cstrehihilo has space for deg+1 doubles;
 *   cstrelohilo has space for deg+1 doubles;
 *   cstrehilolo has space for deg+1 doubles;
 *   cstrelololo has space for deg+1 doubles;
 *   cstimhihihi has space for deg+1 doubles;
 *   cstimlohihi has space for deg+1 doubles;
 *   cstimhilohi has space for deg+1 doubles;
 *   cstimlolohi has space for deg+1 doubles;
 *   cstimhihilo has space for deg+1 doubles;
 *   cstimlohilo has space for deg+1 doubles;
 *   cstimhilolo has space for deg+1 doubles;
 *   cstimlololo has space for deg+1 doubles;
 *   cffrehihihi has space for nbr arrays of deg+1 doubles;
 *   cffrelohihi has space for nbr arrays of deg+1 doubles;
 *   cffrehilohi has space for nbr arrays of deg+1 doubles;
 *   cffrelolohi has space for nbr arrays of deg+1 doubles;
 *   cffrehihilo has space for nbr arrays of deg+1 doubles;
 *   cffrelohilo has space for nbr arrays of deg+1 doubles;
 *   cffrehilolo has space for nbr arrays of deg+1 doubles;
 *   cffrelololo has space for nbr arrays of deg+1 doubles;
 *   cffimhihihi has space for nbr arrays of deg+1 doubles;
 *   cffimlohihi has space for nbr arrays of deg+1 doubles;
 *   cffimhilohi has space for nbr arrays of deg+1 doubles;
 *   cffimlolohi has space for nbr arrays of deg+1 doubles;
 *   cffimhihilo has space for nbr arrays of deg+1 doubles;
 *   cffimlohilo has space for nbr arrays of deg+1 doubles;
 *   cffimhilolo has space for nbr arrays of deg+1 doubles;
 *   cffimlololo has space for nbr arrays of deg+1 doubles;
 *   verbose     if true, then output is written.
 *
 * ON RETURN :
 *   nvr       nvr[k] has the number of variables in monomial k;
 *   idx       idx[k] holds nvr[k] indices to variables in monomial k;
 *   exp       exp[k] holds nvr[k] exponents of variables in monomial k;
 *   inputrehihihi has the highest doubles of dim real series of degree deg;
 *   inputrelohihi has the second highest doubles of dim real input series;
 *   inputrehilohi has the third highest doubles of dim real input series;
 *   inputrelolohi has the fourth highest doubles of dim real input series;
 *   inputrehihilo has the fourth lowest doubles of dim real input series;
 *   inputrelohilo has the third lowest doubles of dim real input series;
 *   inputrehilolo has the second lowest doubles of dim real input series;
 *   inputrelololo has the lowest doubles of dim real series of degree deg;
 *   inputimhihihi has the highest doubles of dim imag series of degree deg;
 *   inputimlohihi has the second highest doubles of dim imag input series;
 *   inputimhilohi has the third highest doubles of dim imag input series;
 *   inputimlolohi has the fourth highest doubles of dim imag input series;
 *   inputimhihilo has the fourth lowest doubles of dim imag input series;
 *   inputimlohilo has the third lowest doubles of dim imag input series;
 *   inputimhilolo has the second lowest doubles of dim imag input series;
 *   inputimlololo has the lowest doubles of dim series of degree deg;
 *   cstrehihihi has the real highest doubles of the constant series;
 *   cstrelohihi has the real second highest doubles of the constant series;
 *   cstrehilohi has the real third highest doubles of the constant series;
 *   cstrelolohi has the real fourth highest doubles of the constant series;
 *   cstrehihilo has the real fourth lowest doubles of the constant series;
 *   cstrelohilo has the real third lowest doubles of the constant series;
 *   cstrehilolo has the real second lowest doubles of the constant series;
 *   cstrelololo has the real lowest doubles of the constant series;
 *   cstimhihihi has the imag highest doubles of the constant series;
 *   cstimlohihi has the imag second highest doubles of the constant series;
 *   cstimhilohi has the imag third highest doubles of the constant series;
 *   cstimlolohi has the imag fourth highest doubles of the constant series;
 *   cstimhihilo has the imag fourth lowest doubles of the constant series;
 *   cstimlohilo has the imag third lowest doubles of the constant series;
 *   cstimhilolo has the imag second lowest doubles of the constant series;
 *   cstimlololo has the imag lowest doubles of the constant series;
 *   cffrehihihi has the highest real doubles of the coefficient
 *               series of monomials;
 *   cffrelohihi has the second real highest doubles of the
 *               coefficient series of monomials;
 *   cffrehilohi has the third real highest doubles of the
 *               coefficient series of monomials;
 *   cffrelolohi has the fourth real highest doubles of the
 *               coefficient series of monomials;
 *   cffrehihilo has the fourth real lowest doubles of the
 *               coefficient series of monomials;
 *   cffrelohilo has the third real lowest doubles of the
 *               coefficient series of monomials;
 *   cffrehilolo has the second real lowest doubles of the
 *               coefficient series of monomials;
 *   cffrelololo has the lowest real doubles of the coefficient
 *               series of monomials.
 *   cffimhihihi has the highest imag doubles of the coefficient
 *               series of monomials;
 *   cffimlohihi has the second imag highest doubles of the
 *               coefficient series of monomials;
 *   cffimhilohi has the third imag highest doubles of the
 *               coefficient series of monomials;
 *   cffimlolohi has the fourth imag highest doubles of the
 *               coefficient series of monomials;
 *   cffimhihilo has the fourth imag lowest doubles of the
 *               coefficient series of monomials;
 *   cffimlohilo has the third imag lowest doubles of the
 *               coefficient series of monomials;
 *   cffimhilolo has the second imag lowest doubles of the
 *               coefficient series of monomials;
 *   cffimlololo has the lowest imag doubles of the coefficient
 *               series of monomials. */

double dbl8_error_sum1
 ( int dim, int deg,
   double **resultshihihi_h, double **resultslohihi_h, 
   double **resultshilohi_h, double **resultslolohi_h,
   double **resultshihilo_h, double **resultslohilo_h, 
   double **resultshilolo_h, double **resultslololo_h,
   double **resultshihihi_d, double **resultslohihi_d,
   double **resultshilohi_d, double **resultslolohi_d,
   double **resultshihilo_d, double **resultslohilo_d,
   double **resultshilolo_d, double **resultslololo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on real data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   resultshihihi_h are the highest doubles computed on the host
 *            without jobs;
 *   resultslohihi_h are the second highest doubles computed on the host
 *            without jobs;
 *   resultshilohi_h are the third highest doubles computed on the host
 *            without jobs;
 *   resultslolohi_h are the fourth highest doubles computed on the host
 *            without jobs;
 *   resultshihilo_h are the fourth lowest doubles computed on the host
 *            without jobs;
 *   resultslohilo_h are the third lowest doubles computed on the host
 *            without jobs;
 *   resultshilolo_h are the second lowest doubles computed on the host
 *            without jobs;
 *   resultslololo_h are the lowest doubles computed on the host
 *            without jobs;
 *   resultshihihi_d are the highest doubles computed on the device;
 *   resultslohihi_d are the second highest doubles computed on the device;
 *   resultshilohi_d are the third highest doubles computed on the device;
 *   resultslolohi_d are the fourth highest doubles computed on the device;
 *   resultshihilo_d are the fourth lowest doubles computed on the device;
 *   resultslohilo_d are the third lowest doubles computed on the device;
 *   resultshilolo_d are the second lowest doubles computed on the device;
 *   resultslololo_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double dbl8_error_sum
 ( int dim, int deg,
   double **results1hihihi_h, double **results1lohihi_h, 
   double **results1hilohi_h, double **results1lolohi_h,
   double **results1hihilo_h, double **results1lohilo_h, 
   double **results1hilolo_h, double **results1lololo_h,
   double **results2hihihi_h, double **results2lohihi_h,
   double **results2hilohi_h, double **results2lolohi_h,
   double **results2hihilo_h, double **results2lohilo_h,
   double **results2hilolo_h, double **results2lololo_h,
   double **resultshihihi_d, double **resultslohihi_d,
   double **resultshilohi_d, double **resultslolohi_d,
   double **resultshihilo_d, double **resultslohilo_d,
   double **resultshilolo_d, double **resultslololo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on real data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   results1hihihi_h are the highest doubles computed on the host
 *            without jobs;
 *   results1lohihi_h are the second highest doubles computed on the host
 *            without jobs;
 *   results1hilohi_h are the third highest doubles computed on the host
 *            without jobs;
 *   results1lolohi_h are the fourth highest doubles computed on the host
 *            without jobs;
 *   results1hihilo_h are the fourth lowest doubles computed on the host
 *            without jobs;
 *   results1lohilo_h are the third lowest doubles computed on the host
 *            without jobs;
 *   results1hilolo_h are the second lowest doubles computed on the host
 *            without jobs;
 *   results1lololo_h are the lowest doubles computed on the host
 *            without jobs;
 *   results2hihihi_h are the highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lohihi_h are the second highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hilohi_h are the third highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lolohi_h are the fourth highest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hihilo_h are the fourth lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lohilo_h are the third lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2hilolo_h are the second lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   results2lololo_h are the lowest doubles computed on the host
 *            with convolution and addition jobs;
 *   resultshihihi_d are the highest doubles computed on the device;
 *   resultslohihi_d are the second highest doubles computed on the device;
 *   resultshilohi_d are the third highest doubles computed on the device;
 *   resultslolohi_d are the fourth highest doubles computed on the device;
 *   resultshihilo_d are the fourth lowest doubles computed on the device;
 *   resultslohilo_d are the third lowest doubles computed on the device;
 *   resultshilolo_d are the second lowest doubles computed on the device;
 *   resultslololo_d are the lowest doubles computed on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double cmplx8_error_sum
 ( int dim, int deg,
   double **resultsrehihihi_h, double **resultsrelohihi_h, 
   double **resultsrehilohi_h, double **resultsrelolohi_h,
   double **resultsrehihilo_h, double **resultsrelohilo_h, 
   double **resultsrehilolo_h, double **resultsrelololo_h,
   double **resultsimhihihi_h, double **resultsimlohihi_h, 
   double **resultsimhilohi_h, double **resultsimlolohi_h,
   double **resultsimhihilo_h, double **resultsimlohilo_h, 
   double **resultsimhilolo_h, double **resultsimlololo_h,
   double **resultsrehihihi_d, double **resultsrelohihi_d,
   double **resultsrehilohi_d, double **resultsrelolohi_d,
   double **resultsrehihilo_d, double **resultsrelohilo_d,
   double **resultsrehilolo_d, double **resultsrelololo_d,
   double **resultsimhihihi_d, double **resultsimlohihi_d,
   double **resultsimhilohi_d, double **resultsimlolohi_d,
   double **resultsimhihilo_d, double **resultsimlohilo_d,
   double **resultsimhilolo_d, double **resultsimlololo_d, bool verbose );
/*
 * DESCRIPTION :
 *   Returns the sum of all errors, comparing results computed on the host
 *   with results computed on the device, on complex data.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   deg      truncation degree of the series;
 *   resultsrehihihi_h are the highest real doubles computed on the host;
 *   resultsrelohihi_h are the 2nd highest real doubles computed on the host;
 *   resultsrehilohi_h are the 3rd highest real doubles computed on the host;
 *   resultsrelolohi_h are the 4th highest real doubles computed on the host;
 *   resultsrehihilo_h are the 4th lowest real doubles computed on the host;
 *   resultsrelohilo_h are the 3rd lowest real doubles computed on the host;
 *   resultsrehilolo_h are the 2nd lowest real doubles computed on the host;
 *   resultsrelololo_h are the lowest real doubles computed on the host;
 *   imsultsimhihihi_h aim the highest imag doubles computed on the host;
 *   imsultsimlohihi_h aim the 2nd highest imag doubles computed on the host;
 *   imsultsimhilohi_h aim the 3rd highest imag doubles computed on the host;
 *   imsultsimlolohi_h aim the 4th highest imag doubles computed on the host;
 *   imsultsimhihilo_h aim the 4th lowest imag doubles computed on the host;
 *   imsultsimlohilo_h aim the 3rd lowest imag doubles computed on the host;
 *   imsultsimhilolo_h aim the 2nd lowest imag doubles computed on the host;
 *   imsultsimlololo_h aim the lowest imag doubles computed on the host;
 *   resultsrehihihi_d are the highest real doubles computed on the device;
 *   resultsrelohihi_d are the 2nd highest real doubles on the device;
 *   resultsrehilohi_d are the 3rd highest real doubles on the device;
 *   resultsrelolohi_d are the 4th highest real doubles on the device;
 *   resultsrehihilo_d are the 4th lowest real doubles on the device;
 *   resultsrelohilo_d are the 3rd lowest real doubles on the device;
 *   resultsrehilolo_d are the 2nd lowest real doubles on the device;
 *   resultsrelololo_d are the lowest real doubles on the device;
 *   imsultsimhihihi_d aim the highest imag doubles computed on the device;
 *   imsultsimlohihi_d aim the 2nd highest imag doubles on the device;
 *   imsultsimhilohi_d aim the 3rd highest imag doubles on the device;
 *   imsultsimlolohi_d aim the 4th highest imag doubles on the device;
 *   imsultsimhihilo_d aim the 4th lowest imag doubles on the device;
 *   imsultsimlohilo_d aim the 3rd lowest imag doubles on the device;
 *   imsultsimhilolo_d aim the 2nd lowest imag doubles on the device;
 *   imsultsimlololo_d aim the lowest imag doubles on the device;
 *   verbose  if true, then all results and intermediate errors are shown. */

double test_dbl8_polynomial
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

double test_cmplx8_polynomial
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
