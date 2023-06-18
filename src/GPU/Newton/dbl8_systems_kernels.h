// The file dbl8_systems_kernels.h specifies functions to define
// memory transfers and kernel launches to evaluate and differentiate
// monomials with common factors in octo double precision.

#ifndef __dbl8_systems_kernels_h__
#define __dbl8_systems_kernels_h__

#include "convolution_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"

void write_dbl8_cnvflops
 ( int dim, int deg, int ctype,
   ConvolutionJobs cnvjobs, double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the kernel and wall time flops for the convolution jobs.
 *
 * ON ENTRY :
 *   dim      number of series;
 *   deg      order of the series, truncation degree;
 *   ctype    0 if on real data, 1 if on complex data;
 *   cnvjobs  defines the convolution jobs;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void write_vectorized8_cnvincflops
 ( int dim, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the kernel and wall time flops for the convolution jobs.
 *
 * ON ENTRY :
 *   dim      number of series;
 *   deg      order of the series, truncation degree;
 *   cnvjobs  defines the convolution jobs;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void dbl8_evaldiffdata_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *
 * ON ENTRY :
 *   datahihihi has the highest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalohihi has the second highest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   datahilohi has the third highest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalolohi has the fourth highest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   datahihilo has the fourth lowest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalohilo has the third lowest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   datahilolo has the second lowest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalololo has the lowest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   outputhihihi has space for the highest doubles of the output;
 *   outputlohihi has space for the second highest doubles of the output;
 *   outputhilohi has space for the third highest doubles of the output;
 *   outputlolohi has space for the fourth highest doubles of the output;
 *   outputhihilo has space for the fourth lowest doubles of the output;
 *   outputlohilo has space for the third lowest doubles of the output;
 *   outputhilolo has space for the second lowest doubles of the output;
 *   outputlololo has space for the lowest doubles of the output;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputhihihi has in outputhihi[i][dim] the highest double value
 *            of the i-th monomial, and
 *            outputhihihi[i][k] has the highest double value
 *            of the k-th derivative of the i-th monomial;
 *   outputlohihi has in outputlohihi[i][dim] the second highest double
 *            value of the i-th monomial, and
 *            outputlohihi[i][k] has the second highest double value
 *            of the k-th derivative of the i-th monomial;
 *   outputhilohi has in outputhilohi[i][dim] the third highest double
 *            value of the i-th monomial, and
 *            outputhilohi[i][k] has the third highest double value
 *            of the k-th derivative of the i-th monomial;
 *   outputlolohi has in outputlolohi[i][dim] the fourth highest double
 *            value of the i-th monomial, and
 *            outputlolohi[i][k] has the fourth highest double value
 *            of the k-th derivative of the i-th monomial;
 *   outputhihilo has in outputhihilo[i][dim] the fourth lowest double
 *            value of the i-th monomial, and
 *            outputhihilo[i][k] has the fourth lowest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputlohilo has in outputlohilo[i][dim] the third lowest double
 *            value of the i-th monomial, and
 *            outputlohilo[i][k] has the third lowest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputhilolo has in outputhilolo[i][dim] the second lowest double
 *            value of the i-th monomial, and
 *            outputhilo[i][k] has the second lowest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputlololo has in outputlololo[i][dim] the lowest double value
 *            of the i-th monomial, and
 *            outputlololo[i][k] has the lowest double value of
 *            the k-th derivative of the i-th monomial. */

void cmplx8_evaldiffdata_to_output
 ( double *datarehihihi, double *datarelohihi,
   double *datarehilohi, double *datarelolohi,
   double *datarehihilo, double *datarelohilo,
   double *datarehilolo, double *datarelololo,
   double *dataimhihihi, double *dataimlohihi,
   double *dataimhilohi, double *dataimlolohi,
   double *dataimhihilo, double *dataimlohilo,
   double *dataimhilolo, double *dataimlololo,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *
 * ON ENTRY :
 *   datarehihihi are the highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelohihi are the 2nd highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarehilohi are the 3rd highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelolohi are the 4th highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarehihilo are the 4th lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelohilo are the 3rd lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarehilolo are the 2nd lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelololo are the lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   dataimhihihi are the highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlohihi are the 2nd highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimhilohi are the 3rd highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlolohi are the 4th highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimhilhio are the 4th lowest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlohilo are the 3rd lowest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimhilolo are the 2nd lowest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlololo are the lowest doubles of the imaginary parts of the 
 *            input, computed forward, backward, and cross products;
 *   outputrehihihi has space for the value and all derivatives;
 *   outputrelohihi has space for the value and all derivatives;
 *   outputrehilohi has space for the value and all derivatives;
 *   outputrelolohi has space for the value and all derivatives;
 *   outputrehihilo has space for the value and all derivatives;
 *   outputrelohilo has space for the value and all derivatives;
 *   outputrehilolo has space for the value and all derivatives;
 *   outputrelololo has space for the value and all derivatives;
 *   outputimhihihi has space for the value and all derivatives;
 *   outputimlohihi has space for the value and all derivatives;
 *   outputimhilohi has space for the value and all derivatives;
 *   outputimlolohi has space for the value and all derivatives;
 *   outputimhihilo has space for the value and all derivatives;
 *   outputimlohilo has space for the value and all derivatives;
 *   outputimhilolo has space for the value and all derivatives;
 *   outputimlololo has space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihihi, in outputrehihihi[i][dim] are the highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrehihihi[i][k] has the highest doubles of the real parts
 *            of the value of the k-th derivative of the i-th monomial;
 *   outputrelohihi, in outputrelohi[i][dim] are the 2nd highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelohihi[i][k] has the 2nd highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilohi, in outputrelohi[i][dim] are the 3rd highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrehilohi[i][k] has the 3rd highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelolohi, in outputrelolohi[i][dim] are the 4th highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelolohi[i][k] has the 4th highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehihilo, in outputrehihilo[i][dim] are the 4th lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehihilo[i][k] has the 4th lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelohilo, in outputrelohilo[i][dim] are the 3rd lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrelohilo[i][k] has the 3rd lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilolo, in outputrehilolo[i][dim] are the 2nd lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehilolo[i][k] has the 2nd lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelololo, in outputrelololo[i][dim] are the lowest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelolo[i][k] has the lowest doubles of the real parts
 *            of the value of the k-th derivative of the i-th monomial;
 *   outputimhihihi, in outputimhihihi[i][dim] are the highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhihihi[i][k] has the highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohihi, in outputrelohihi[i][dim] are the 2nd highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohihi[i][k] has the 2nd highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilohi, in outputrehilohi[i][dim] are the 3nd highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhilohi[i][k] has the 3rd highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlolohi, in outputrelolohi[i][dim] are the 4th highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolohi[i][k] has the 4th highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhihilo, in outputrehihilo[i][dim] are the 4th lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolohi[i][k] has the 4th lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohilo, in outputrelohilo[i][dim] are the 3rd lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohilo[i][k] has the 3rd lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilolo, in outputrehilolo[i][dim] are the 2nd lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhilohi[i][k] has the 2nd lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlololo, in outputrelololo[i][dim] are the lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlololo[i][k] has the lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial. */

void cmplx8vectorized_evaldiffdata_to_output
 ( double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totalcff, int offsetri, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device
 *   using vectorized arithmetic to the output.
 *
 * ON ENTRY :
 *   datarihihihi are the highest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarilohihi are the 2nd highest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarihilohi are the 3rd highest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarilolohi are the 4th highest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarihihilo are the 4th lowest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarilohilo are the 3rd lowest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarihilolo are the 2nd lowest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   datarilololo are the lowest doubles of the input, 
 *            computed forward, backward, and cross products;
 *   outputrehihihi has space for the value and all derivatives;
 *   outputrelohihi has space for the value and all derivatives;
 *   outputrehilohi has space for the value and all derivatives;
 *   outputrelolohi has space for the value and all derivatives;
 *   outputrehihilo has space for the value and all derivatives;
 *   outputrelohilo has space for the value and all derivatives;
 *   outputrehilolo has space for the value and all derivatives;
 *   outputrelololo has space for the value and all derivatives;
 *   outputimhihihi has space for the value and all derivatives;
 *   outputimlohihi has space for the value and all derivatives;
 *   outputimhilohi has space for the value and all derivatives;
 *   outputimlolohi has space for the value and all derivatives;
 *   outputimhihilo has space for the value and all derivatives;
 *   outputimlohilo has space for the value and all derivatives;
 *   outputimhilolo has space for the value and all derivatives;
 *   outputimlololo has space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   totalcff is the total number of coefficients without vectorization;
 *   offsetri is the size of the second operand;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihihi, in outputrehihihi[i][dim] are the highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrehihihi[i][k] has the highest doubles of the real parts
 *            of the value of the k-th derivative of the i-th monomial;
 *   outputrelohihi, in outputrelohi[i][dim] are the 2nd highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelohihi[i][k] has the 2nd highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilohi, in outputrelohi[i][dim] are the 3rd highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrehilohi[i][k] has the 3rd highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelolohi, in outputrelolohi[i][dim] are the 4th highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelolohi[i][k] has the 4th highest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehihilo, in outputrehihilo[i][dim] are the 4th lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehihilo[i][k] has the 4th lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelohilo, in outputrelohilo[i][dim] are the 3rd lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrelohilo[i][k] has the 3rd lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilolo, in outputrehilolo[i][dim] are the 2nd lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehilolo[i][k] has the 2nd lowest doubles of the real
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrelololo, in outputrelololo[i][dim] are the lowest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelolo[i][k] has the lowest doubles of the real parts
 *            of the value of the k-th derivative of the i-th monomial;
 *   outputimhihihi, in outputimhihihi[i][dim] are the highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhihihi[i][k] has the highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohihi, in outputrelohihi[i][dim] are the 2nd highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohihi[i][k] has the 2nd highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilohi, in outputrehilohi[i][dim] are the 3nd highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhilohi[i][k] has the 3rd highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlolohi, in outputrelolohi[i][dim] are the 4th highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolohi[i][k] has the 4th highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhihilo, in outputrehihilo[i][dim] are the 4th lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolohi[i][k] has the 4th lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohilo, in outputrelohilo[i][dim] are the 3rd lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohilo[i][k] has the 3rd lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilolo, in outputrehilolo[i][dim] are the 2nd lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhilohi[i][k] has the 2nd lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlololo, in outputrelololo[i][dim] are the lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlololo[i][k] has the lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial. */

void GPU_dbl8_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   on real data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cffhihihi, cffhihihi[k] has deg+1 doubles for the highest doubles
 *            of the coefficient series of monomial k;
 *   cfflohihi, cfflohihi[k] has deg+1 doubles for the second highest doubles
 *            of the coefficient series of monomial k;
 *   cffhilohi, cffhilohi[k] has deg+1 doubles for the third highest doubles
 *            of the coefficient series of monomial k;
 *   cfflolohi, cfflolohi[k] has deg+1 doubles for the fourth highest doubles
 *            of the coefficient series of monomial k;
 *   cffhihilo, cffhihilo[k] has deg+1 doubles for the fourth lowest doubles
 *            of the coefficient series of monomial k;
 *   cfflohilo, cfflohilo[k] has deg+1 doubles for the third lowest doubles
 *            of the coefficient series of monomial k;
 *   cffhilolo, cffhilolo[k] has deg+1 doubles for the second lowest doubles
 *            of the coefficient series of monomial k;
 *   cfflololo, cfflololo[k] has deg+1 doubles for the lowest doubles
 *            of the coefficient series of monomial k;
 *   inputhihihi has the highest doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputlohihi has the second highest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputhilohi has the third highest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputlolohi has the fourth highest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputhihilo has the fourth lowest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputlohilo has the third lowest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputhilolo has the second lowest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   inputlololo has the lowest doubles of the coefficients
 *            of the series for all variables in the polynomial;
 *   outputhihihi has space allocated for the highest double output;
 *   outputlohihi has space allocated for the second highest double output;
 *   outputhilohi has space allocated for the third highest double output;
 *   outputlolohi has space allocated for the fourth highest double output;
 *   outputhihilo has space allocated for the fourth lowest double output;
 *   outputlohilo has space allocated for the third lowest double output;
 *   outputhilolo has space allocated for the second lowest double output;
 *   outputlololo has space allocated for the lowest double output;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputhihihi has the highest doubles of the output,
 *            outputhihihi[k], for k from 0 to dim-1, has the highest double
 *            of the derivative with respect to the variable k;
 *            outputhihihi[dim] has the highest double value;
 *   outputlohihi has the second highest doubles of the output,
 *            outputlohihi[k], for k from 0 to dim-1, has the second highest
 *            double of the derivative with respect to the variable k;
 *            outputlohihi[dim] has the second highest double value;
 *   outputhilohi has the second lowest doubles of the output,
 *            outputhilohi[k], for k from 0 to dim-1, has the third highest
 *            double of the derivative with respect to the variable k;
 *            outputhilohi[dim] has the third highest double value;
 *   outputlolohi has the fourth highest doubles of the output,
 *            outputlolohi[k], for k from 0 to dim-1, has the fourth highest
 *            double of the derivative with respect to the variable k;
 *            outputlolohi[dim] has the fourth highest double value;
 *   outputhihilo has the fourth lowest doubles of the output,
 *            outputhihilo[k], for k from 0 to dim-1, has the fourth lowest
 *            double of the derivative with respect to the variable k;
 *            outputhihilo[dim] has the fourth lowest double value;
 *   outputlohilo has the third lowest doubles of the output,
 *            outputlohilo[k], for k from 0 to dim-1, has the third lowest
 *            double of the derivative with respect to the variable k;
 *            outputlohilo[dim] has the third lowest double value;
 *   outputhilolo has the second lowest doubles of the output,
 *            outputhilolo[k], for k from 0 to dim-1, has the second lowest
 *            double of the derivative with respect to the variable k;
 *            outputhilolo[dim] has the second lowest double value;
 *   outputlololo has the lowest doubles of the output,
 *            outputlololo[k], for k from 0 to dim-1, has the lowest double
 *            of the derivative with respect to the variable k;
 *            outputlololo[dim] has the lowest double value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx8_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   on complex data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cffrehihihi, cffrehihihi[k] has deg+1 highest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffrelohihi, cffrelohihi[k] has deg+1 second highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilohi, cffrehilohi[k] has deg+1 third highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelolohi, cffrelolohi[k] has deg+1 fourth highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehihilo, cffrehihilo[k] has deg+1 fourth lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelohilo, cffrelohilo[k] has deg+1 third lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilolo, cffrehilolo[k] has deg+1 second lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelololo, cffrelololo[k] has deg+1 lowest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffimhihihi, cffimhihihi[k] has deg+1 highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohihi, cffimlohihi[k] has deg+1 second highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilohi, cffimhilohi[k] has deg+1 third highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlolohi, cffimlolohi[k] has deg+1 fourth highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhihilo, cffimihihilo[k] has deg+1 fourth lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohilo, cffimilohilo[k] has deg+1 third lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilolo, cffimihilolo[k] has deg+1 second lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlololo, cffimlololo[k] has deg+1 lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   inputrehihihi are the highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohihi are the second highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilohi are the third highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelolohi are the fourth highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehihilo are the fourth lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohilo are the third lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilolo are the second lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelololo are the lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihihi are the highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohihi are the second highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilohi are the third highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlolohi are the fourth lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihilo are the fourth lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohilo are the third lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilolo are the second lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlololo are the lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   outputrehihihi has space for the highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohihi has space for the second highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilohi has space for the third highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelolohi has space for the fourth highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehihilo has space for the fourth lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohilo has space for the third lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilolo has space for the second lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelololo has space for the lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputimhihihi has space for the highest doubles of the imaginary parts
 *            of value and all derivatives;
 *   outputimlohihi has space for the second highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilohi has space for the third highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlolohi has space for the fourth highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhihilo has space for the fourth lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlohilo has space for the third lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilolo has space for the second lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlololo has space for the lowest doubles of the imaginary parts
 *            of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihihi has the highest doubles of the real parts of derivatives
 *            and the value, outputrehihihi[k], for k from 0 to dim-1,
 *            has the highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihihi[dim] has the
 *            highest doubles of the real part of the value;
 *   outputrelohihi has the 2nd highest doubles of the real parts of derivatives
 *            and the value, outputrelohihi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohihi[dim] has the
 *            second highest doubles of the real part of the value;
 *   outputrehilohi has the 3rd highest doubles of the real parts of derivatives
 *            and the value, outputrehilohi[k], for k from 0 to dim-1, has the
 *            third highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilohi[dim] has the
 *            third highest doubles of the real part of the value;
 *   outputrelolohi has the 4th highest doubles of the real parts of derivatives
 *            and the value, outputrelolohi[k], for k from 0 to dim-1, has the
 *            fourth highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelolohi[dim] has the
 *            fourth highest doubles of the real part of the value;
 *   outputrehihilo has the 4th lowest doubles of the real parts of derivatives
 *            and the value, outputrehihilo[k], for k from 0 to dim-1, has the
 *            fourth lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihilo[dim] has the
 *            fourth lowest doubles of the real part of the value;
 *   outputrelohilo has the 3rd lowest doubles of the real parts of derivatives
 *            and the value, outputrelohilo[k], for k from 0 to dim-1, has the
 *            third lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohilo[dim] has the
 *            third lowest doubles of the real part of the value;
 *   outputrehilolo has the 2nd lowest doubles of the real parts of derivatives
 *            and the value, outputrehilolo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilolo[dim] has the
 *            second lowest doubles of the real part of the value;
 *   outputrelololo has the lowest doubles of the real parts of derivatives
 *            and the value, outputrelololo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelololo[dim] has the
 *            lowest doubles of the real part of the value;
 *   outputimhihihi has the highest doubles of the imag parts of derivatives
 *            and the value, outputimhihihi[k], for k from 0 to dim-1, has the
 *            highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihihi[dim] has the
 *            highest doubles of the imaginary part of the value;
 *   outputimlohihi has the 2nd highest doubles of the imag parts of derivatives
 *            and the value, outputimlohihi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohihi[dim] has the
 *            second highest doubles of the imaginary part of the value;
 *   outputimhilohi has the 3rd highest doubles of the imag parts of derivatives
 *            and the value, outputimhilohi[k], for k from 0 to dim-1, has the
 *            third highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilohi[dim] has the
 *            third highest doubles of the imaginary part of the value;
 *   outputimlolohi has the 4th highest doubles of the imag parts of derivatives
 *            and the value, outputimlolohi[k], for k from 0 to dim-1, has the
 *            fourth highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlolohi[dim] has the
 *            fourth highest doubles of the imaginary part of the value;
 *   outputimhihilo has the 4th lowest doubles of the imag parts of derivatives
 *            and the value, outputimhihilo[k], for k from 0 to dim-1, has the
 *            fourth lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihilo[dim] has the
 *            fourth lowest doubles of the imaginary part of the value;
 *   outputimlohilo has the 3rd lowest doubles of the imag parts of derivatives
 *            and the value, outputimlohilo[k], for k from 0 to dim-1, has the
 *            third lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohilo[dim] has the
 *            third lowest doubles of the imaginary part of the value;
 *   outputimhilolo has the 2nd lowest doubles of the imag parts of derivatives
 *            and the value, outputimhilolo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilolo[dim] has the
 *            second lowest doubles of the imaginary part of the value;
 *   outputimlololo has the lowest doubles of the imaginary parts of derivatives
 *            and the value, outputimlololo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlololo[dim] has the
 *            lowest doubles of the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx8vectorized_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   on complex data using complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cffrehihihi, cffrehihihi[k] has deg+1 highest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffrelohihi, cffrelohihi[k] has deg+1 second highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilohi, cffrehilohi[k] has deg+1 third highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelolohi, cffrelolohi[k] has deg+1 fourth highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehihilo, cffrehihilo[k] has deg+1 fourth lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelohilo, cffrelohilo[k] has deg+1 third lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilolo, cffrehilolo[k] has deg+1 second lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelololo, cffrelololo[k] has deg+1 lowest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffimhihihi, cffimhihihi[k] has deg+1 highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohihi, cffimlohihi[k] has deg+1 second highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilohi, cffimhilohi[k] has deg+1 third highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlolohi, cffimlolohi[k] has deg+1 fourth highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhihilo, cffimihihilo[k] has deg+1 fourth lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohilo, cffimilohilo[k] has deg+1 third lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilolo, cffimihilolo[k] has deg+1 second lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlololo, cffimlololo[k] has deg+1 lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   inputrehihihi are the highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohihi are the second highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilohi are the third highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelolohi are the fourth highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehihilo are the fourth lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohilo are the third lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilolo are the second lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelololo are the lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihihi are the highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohihi are the second highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilohi are the third highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlolohi are the fourth lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihilo are the fourth lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohilo are the third lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilolo are the second lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlololo are the lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   outputrehihihi has space for the highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohihi has space for the second highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilohi has space for the third highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelolohi has space for the fourth highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehihilo has space for the fourth lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohilo has space for the third lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilolo has space for the second lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelololo has space for the lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputimhihihi has space for the highest doubles of the imaginary parts
 *            of value and all derivatives;
 *   outputimlohihi has space for the second highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilohi has space for the third highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlolohi has space for the fourth highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhihilo has space for the fourth lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlohilo has space for the third lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilolo has space for the second lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlololo has space for the lowest doubles of the imaginary parts
 *            of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihihi has the highest doubles of the real parts of derivatives
 *            and the value, outputrehihihi[k], for k from 0 to dim-1,
 *            has the highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihihi[dim] has the
 *            highest doubles of the real part of the value;
 *   outputrelohihi has the 2nd highest doubles of the real parts of derivatives
 *            and the value, outputrelohihi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohihi[dim] has the
 *            second highest doubles of the real part of the value;
 *   outputrehilohi has the 3rd highest doubles of the real parts of derivatives
 *            and the value, outputrehilohi[k], for k from 0 to dim-1, has the
 *            third highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilohi[dim] has the
 *            third highest doubles of the real part of the value;
 *   outputrelolohi has the 4th highest doubles of the real parts of derivatives
 *            and the value, outputrelolohi[k], for k from 0 to dim-1, has the
 *            fourth highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelolohi[dim] has the
 *            fourth highest doubles of the real part of the value;
 *   outputrehihilo has the 4th lowest doubles of the real parts of derivatives
 *            and the value, outputrehihilo[k], for k from 0 to dim-1, has the
 *            fourth lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihilo[dim] has the
 *            fourth lowest doubles of the real part of the value;
 *   outputrelohilo has the 3rd lowest doubles of the real parts of derivatives
 *            and the value, outputrelohilo[k], for k from 0 to dim-1, has the
 *            third lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohilo[dim] has the
 *            third lowest doubles of the real part of the value;
 *   outputrehilolo has the 2nd lowest doubles of the real parts of derivatives
 *            and the value, outputrehilolo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilolo[dim] has the
 *            second lowest doubles of the real part of the value;
 *   outputrelololo has the lowest doubles of the real parts of derivatives
 *            and the value, outputrelololo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelololo[dim] has the
 *            lowest doubles of the real part of the value;
 *   outputimhihihi has the highest doubles of the imag parts of derivatives
 *            and the value, outputimhihihi[k], for k from 0 to dim-1, has the
 *            highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihihi[dim] has the
 *            highest doubles of the imaginary part of the value;
 *   outputimlohihi has the 2nd highest doubles of the imag parts of derivatives
 *            and the value, outputimlohihi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohihi[dim] has the
 *            second highest doubles of the imaginary part of the value;
 *   outputimhilohi has the 3rd highest doubles of the imag parts of derivatives
 *            and the value, outputimhilohi[k], for k from 0 to dim-1, has the
 *            third highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilohi[dim] has the
 *            third highest doubles of the imaginary part of the value;
 *   outputimlolohi has the 4th highest doubles of the imag parts of derivatives
 *            and the value, outputimlolohi[k], for k from 0 to dim-1, has the
 *            fourth highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlolohi[dim] has the
 *            fourth highest doubles of the imaginary part of the value;
 *   outputimhihilo has the 4th lowest doubles of the imag parts of derivatives
 *            and the value, outputimhihilo[k], for k from 0 to dim-1, has the
 *            fourth lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihilo[dim] has the
 *            fourth lowest doubles of the imaginary part of the value;
 *   outputimlohilo has the 3rd lowest doubles of the imag parts of derivatives
 *            and the value, outputimlohilo[k], for k from 0 to dim-1, has the
 *            third lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohilo[dim] has the
 *            third lowest doubles of the imaginary part of the value;
 *   outputimhilolo has the 2nd lowest doubles of the imag parts of derivatives
 *            and the value, outputimhilolo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilolo[dim] has the
 *            second lowest doubles of the imaginary part of the value;
 *   outputimlololo has the lowest doubles of the imaginary parts of derivatives
 *            and the value, outputimlololo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlololo[dim] has the
 *            lowest doubles of the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_dbl8_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim      number of monomials;
 *   deg      degree of the power series;
 *   szt      size of each block of threads;
 *   nbt      number of thread blocks;
 *   nvr      nvr[i] is the number of variables in the i-th monomial;
 *   idx      idx[i] are the indices of the variables in monomial i;
 *   exp      exp[i] are the exponents of the variables in monomial i;
 *   nbrfac   nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac   expfac[i] are the exponents in the i-th polynomial
 *            that are larger than one, minus one in the factor,
 *            if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffhihihi has the highest doubles of the coefficients of the monomials;
 *   cfflohihi has the second highest doubles of the coefficients;
 *   cffhilohi has the third highest doubles of the coefficients;
 *   cfflolohi has the fourth highest doubles of the coefficients;
 *   cffhihilo has the fourth lowest doubles of the coefficients;
 *   cfflohilo has the third lowest doubles of the coefficients;
 *   cffhilolo has the second lowest doubles of the coefficients;
 *   cfflololo has the lowest doubles of the coefficients of the monomials;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   inputhihi has the highest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputlohi has the second highest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputhilo has the second lowest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputlolo has the lowest doubles of the coefficients of the input,
 *            for dim variables;
 *   outputhihihi has space for the highest doubles of the output;
 *   outputlohihi has space for the second highest doubles of the output;
 *   outputhilohi has space for the third highest doubles of the output;
 *   outputlolohi has space for the fourth highest doubles of the output;
 *   outputhihilo has space for the fourth lowest doubles of the output;
 *   outputlohilo has space for the third lowest doubles of the output;
 *   outputhilolo has space for the second lowest doubles of the output;
 *   outputlololo has space for the lowest doubles of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffhihihi has the highest doubles of the common factors;
 *   cfflohihi has the second highest doubles of the common factors;
 *   cffhilohi has the third highest doubles of the common factors;
 *   cfflolohi has the fourth highest doubles of the common factors;
 *   cffhihilo has the fourth lowest doubles of the common factors;
 *   cfflohilo has the third lowest doubles of the common factors;
 *   cffhilolo has the second lowest doubles of the common factors;
 *   cfflololo has the lowest doubles of the evaluated common factors;
 *   outputhihihi has the highest doubles of the output,
 *            outputhihihi[k], for k from 0 to dim-1, has the highest double
 *            of the derivative with respect to the variable k;
 *            outputhihihi[dim] has the highest double value;
 *   outputlohihi has the second highest doubles of the output,
 *            outputlohihi[k], for k from 0 to dim-1, has the second highest
 *            double of the derivative with respect to the variable k;
 *            outputlohihi[dim] has the second highest double value;
 *   outputhilohi has the second lowest doubles of the output,
 *            outputhilohi[k], for k from 0 to dim-1, has the third highest
 *            double of the derivative with respect to the variable k;
 *            outputhilohi[dim] has the third highest double value;
 *   outputlolohi has the fourth highest doubles of the output,
 *            outputlolohi[k], for k from 0 to dim-1, has the fourth highest
 *            double of the derivative with respect to the variable k;
 *            outputlolohi[dim] has the fourth highest double value;
 *   outputhihilo has the fourth lowest doubles of the output,
 *            outputhihilo[k], for k from 0 to dim-1, has the fourth lowest
 *            double of the derivative with respect to the variable k;
 *            outputhihilo[dim] has the fourth lowest double value;
 *   outputlohilo has the third lowest doubles of the output,
 *            outputlohilo[k], for k from 0 to dim-1, has the third lowest
 *            double of the derivative with respect to the variable k;
 *            outputlohilo[dim] has the third lowest double value;
 *   outputhilolo has the second lowest doubles of the output,
 *            outputhilolo[k], for k from 0 to dim-1, has the second lowest
 *            double of the derivative with respect to the variable k;
 *            outputhilolo[dim] has the second lowest double value;
 *   outputlololo has the lowest doubles of the output,
 *            outputlololo[k], for k from 0 to dim-1, has the lowest double
 *            of the derivative with respect to the variable k;
 *            outputlololo[dim] has the lowest double value;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx8_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double *accrehihihi, double *accrelohihi,
   double *accrehilohi, double *accrelolohi,
   double *accrehihilo, double *accrelohilo,
   double *accrehilolo, double *accrelololo,
   double *accimhihihi, double *accimlohihi,
   double *accimhilohi, double *accimlolohi,
   double *accimhihilo, double *accimlohilo,
   double *accimhilolo, double *accimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi, 
   double **inputimhilohi, double **inputimlolohi, 
   double **inputimhihilo, double **inputimlohilo, 
   double **inputimhilolo, double **inputimlololo, 
   double ***outputrehihihi, double ***outputrelohihi, 
   double ***outputrehilohi, double ***outputrelolohi, 
   double ***outputrehihilo, double ***outputrelohilo, 
   double ***outputrehilolo, double ***outputrelololo, 
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffrehihihi are the highest doubles of the real parts
 *             of the coefficients;
 *   cffrelohihi are the second highest doubles of the real parts
 *             of the coefficients;
 *   cffrehilohi are the third highest doubles of the real parts
 *             of the coefficients;
 *   cffrelolohi are the fourth highest doubles of the real parts
 *             of the coefficients;
 *   cffrehihilo are the fourth lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelohilo are the third lowest doubles of the real parts
 *             of the coefficients;
 *   cffrehilolo are the second lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelololo are the lowest doubles of the real parts of the coefficients;
 *   cffimhihihi are the highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohihi are the second highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilohi are the third highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlolohi are the fourth highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhihilo are the fourth lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohilo are the third lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilolo are the second lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlololo are the lowest doubles of the imaginary parts
 *             of the coefficients;
 *   accrehihihi has space to accumulate one series of degree deg;
 *   accrelohihi has space to accumulate one series of degree deg;
 *   accrehilohi has space to accumulate one series of degree deg;
 *   accrelolohi has space to accumulate one series of degree deg;
 *   accrehihilo has space to accumulate one series of degree deg;
 *   accrelohilo has space to accumulate one series of degree deg;
 *   accrehilolo has space to accumulate one series of degree deg;
 *   accrelololo has space to accumulate one series of degree deg;
 *   accimhihihi has space to accumulate one  series of degree deg;
 *   accimlohihi has space to accumulate one  series of degree deg;
 *   accimhilohi has space to accumulate one  series of degree deg;
 *   accimlolohi has space to accumulate one  series of degree deg;
 *   accimhihilo has space to accumulate one  series of degree deg;
 *   accimlohilo has space to accumulate one  series of degree deg;
 *   accimhilolo has space to accumulate one  series of degree deg;
 *   accimlololo has space to accumulate one  series of degree deg;
 *   inputrehihihi has the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohihi has the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilohi has the third highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohihi has the fourth highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehihilo has the fourth lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohilo has the third lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilolo has the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelololo has the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihihi has the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohihi has the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilohi has the third highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolohi has the fourth highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihilo has the fourth lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohilo has the third lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilolo has the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlololo has the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihihi has space for the highest doubles of the real parts
 *             of the output;
 *   outputrelohihi has space for the second highest doubles of the real parts
 *             of the output;
 *   outputrehilohi has space for the third highest doubles of the real parts
 *             of the output;
 *   outputrelolohi has space for the fourth highest doubles of the real parts
 *             of the output;
 *   outputrehihilo has space for the fourth lowest doubles of the real parts
 *             of the output;
 *   outputrelohilo has space for the third lowest doubles of the real parts
 *             of the output;
 *   outputrehilolo has space for the second lowest doubles of the real parts
 *             of the output;
 *   outputrelololo has space for the lowest doubles of the real parts
 *             of the output;
 *   outputimhihihi has space for the highest doubles of the imaginary
 *             parts of the output;
 *   outputimlohihi has space for the second highest doubles of the imaginary
 *             parts of the output;
 *   outputimhilohi has space for the third highest doubles of the imaginary
 *             parts of the output;
 *   outputimlolohi has space for the fourth highest doubles of the imaginary
 *             parts of the output;
 *   outputimhihilo has space for the fourth lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlohilo has space for the third lowest doubles of the imaginary
 *             parts of the output;
 *   outputimhilolo has space for the second lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlololo has space for the lowest doubles of the imaginary
 *             parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffrehihihi has the highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohihi has the second highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilohi has the third highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelolohi has the fourth highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehihilo has the fourth lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohilo has the third lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilolo has the second lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelololo has the lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffimhihihi has the highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohihi has the second highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilohi has the third highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlolohi has the fourth highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhihilo has the fourth lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohilo has the third lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilolo has the second lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlololo has the lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   outputrehihihi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrehihihi[i][dim] has the
 *             highest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihihi[i][idx[k]] has the highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohihi has the second highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelohihi[i][dim]
 *             has the second highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohihi[i][idx[k]] has the second highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilohi has the third highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilohi[i][dim]
 *             has the third highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilohi[i][idx[k]] has the third highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelolohi has the fourth highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelolohi[i][dim]
 *             has the fourth highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelolohi[i][idx[k]] has the fourth highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehihilo has the fourth lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehihilo[i][dim]
 *             has the fourth lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihilo[i][idx[k]] has the fourth lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohilo has the third lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelohilo[i][dim]
 *             has the third lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohilo[i][idx[k]] has the third lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilolo has the second lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilolo[i][dim]
 *             has the second lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilolo[i][idx[k]] has the second lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelololo has the lowest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelololo[i][dim] has the
 *             lowest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelololo[i][idx[k]] has the lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputimhihihi has the highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhihihi[i][dim]
 *             has the highest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimhihi[i][idx[k]] has the highest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohihi has the second highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlohihi[i][dim]
 *             has the second highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohihi[i][idx[k]] has the second highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilohi has the third highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhilohi[i][dim]
 *             has the third highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilohi[i][idx[k]] has the third highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlolohi has the fourth highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlolohi[i][dim]
 *             has the fourth highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlolohi[i][idx[k]] has the fourth highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhihilo has the fourth lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhihilo[i][dim]
 *             has the fourth lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhihilo[i][idx[k]] has the fourth lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohilo has the third lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlohilo[i][dim]
 *             has the third lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohilo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilolo has the second lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhilolo[i][dim]
 *             has the second lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilolo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlololo has the lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlololo[i][dim]
 *             has the lowest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimlololo[i][idx[k]] has the lowest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */


void GPU_cmplx8vectorized_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double *accrehihihi, double *accrelohihi,
   double *accrehilohi, double *accrelolohi,
   double *accrehihilo, double *accrelohilo,
   double *accrehilolo, double *accrelololo,
   double *accimhihihi, double *accimlohihi,
   double *accimhilohi, double *accimlolohi,
   double *accimhihilo, double *accimlohilo,
   double *accimhilolo, double *accimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi, 
   double **inputimhilohi, double **inputimlolohi, 
   double **inputimhihilo, double **inputimlohilo, 
   double **inputimhilolo, double **inputimlololo, 
   double ***outputrehihihi, double ***outputrelohihi, 
   double ***outputrehilohi, double ***outputrelolohi, 
   double ***outputrehihilo, double ***outputrelohilo, 
   double ***outputrehilolo, double ***outputrelololo, 
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series using vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffrehihihi are the highest doubles of the real parts
 *             of the coefficients;
 *   cffrelohihi are the second highest doubles of the real parts
 *             of the coefficients;
 *   cffrehilohi are the third highest doubles of the real parts
 *             of the coefficients;
 *   cffrelolohi are the fourth highest doubles of the real parts
 *             of the coefficients;
 *   cffrehihilo are the fourth lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelohilo are the third lowest doubles of the real parts
 *             of the coefficients;
 *   cffrehilolo are the second lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelololo are the lowest doubles of the real parts of the coefficients;
 *   cffimhihihi are the highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohihi are the second highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilohi are the third highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlolohi are the fourth highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhihilo are the fourth lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohilo are the third lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilolo are the second lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlololo are the lowest doubles of the imaginary parts
 *             of the coefficients;
 *   accrehihihi has space to accumulate one series of degree deg;
 *   accrelohihi has space to accumulate one series of degree deg;
 *   accrehilohi has space to accumulate one series of degree deg;
 *   accrelolohi has space to accumulate one series of degree deg;
 *   accrehihilo has space to accumulate one series of degree deg;
 *   accrelohilo has space to accumulate one series of degree deg;
 *   accrehilolo has space to accumulate one series of degree deg;
 *   accrelololo has space to accumulate one series of degree deg;
 *   accimhihihi has space to accumulate one  series of degree deg;
 *   accimlohihi has space to accumulate one  series of degree deg;
 *   accimhilohi has space to accumulate one  series of degree deg;
 *   accimlolohi has space to accumulate one  series of degree deg;
 *   accimhihilo has space to accumulate one  series of degree deg;
 *   accimlohilo has space to accumulate one  series of degree deg;
 *   accimhilolo has space to accumulate one  series of degree deg;
 *   accimlololo has space to accumulate one  series of degree deg;
 *   inputrehihihi has the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohihi has the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilohi has the third highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohihi has the fourth highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehihilo has the fourth lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohilo has the third lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilolo has the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelololo has the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihihi has the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohihi has the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilohi has the third highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolohi has the fourth highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihilo has the fourth lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohilo has the third lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilolo has the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlololo has the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihihi has space for the highest doubles of the real parts
 *             of the output;
 *   outputrelohihi has space for the second highest doubles of the real parts
 *             of the output;
 *   outputrehilohi has space for the third highest doubles of the real parts
 *             of the output;
 *   outputrelolohi has space for the fourth highest doubles of the real parts
 *             of the output;
 *   outputrehihilo has space for the fourth lowest doubles of the real parts
 *             of the output;
 *   outputrelohilo has space for the third lowest doubles of the real parts
 *             of the output;
 *   outputrehilolo has space for the second lowest doubles of the real parts
 *             of the output;
 *   outputrelololo has space for the lowest doubles of the real parts
 *             of the output;
 *   outputimhihihi has space for the highest doubles of the imaginary
 *             parts of the output;
 *   outputimlohihi has space for the second highest doubles of the imaginary
 *             parts of the output;
 *   outputimhilohi has space for the third highest doubles of the imaginary
 *             parts of the output;
 *   outputimlolohi has space for the fourth highest doubles of the imaginary
 *             parts of the output;
 *   outputimhihilo has space for the fourth lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlohilo has space for the third lowest doubles of the imaginary
 *             parts of the output;
 *   outputimhilolo has space for the second lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlololo has space for the lowest doubles of the imaginary
 *             parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffrehihihi has the highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohihi has the second highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilohi has the third highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelolohi has the fourth highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehihilo has the fourth lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohilo has the third lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilolo has the second lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelololo has the lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffimhihihi has the highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohihi has the second highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilohi has the third highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlolohi has the fourth highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhihilo has the fourth lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohilo has the third lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilolo has the second lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlololo has the lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   outputrehihihi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrehihihi[i][dim] has the
 *             highest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihihi[i][idx[k]] has the highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohihi has the second highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelohihi[i][dim]
 *             has the second highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohihi[i][idx[k]] has the second highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilohi has the third highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilohi[i][dim]
 *             has the third highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilohi[i][idx[k]] has the third highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelolohi has the fourth highest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelolohi[i][dim]
 *             has the fourth highest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelolohi[i][idx[k]] has the fourth highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehihilo has the fourth lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehihilo[i][dim]
 *             has the fourth lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihilo[i][idx[k]] has the fourth lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohilo has the third lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrelohilo[i][dim]
 *             has the third lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohilo[i][idx[k]] has the third lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilolo has the second lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilolo[i][dim]
 *             has the second lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilolo[i][idx[k]] has the second lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelololo has the lowest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelololo[i][dim] has the
 *             lowest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelololo[i][idx[k]] has the lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputimhihihi has the highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhihihi[i][dim]
 *             has the highest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimhihi[i][idx[k]] has the highest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohihi has the second highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlohihi[i][dim]
 *             has the second highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohihi[i][idx[k]] has the second highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilohi has the third highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhilohi[i][dim]
 *             has the third highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilohi[i][idx[k]] has the third highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlolohi has the fourth highest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlolohi[i][dim]
 *             has the fourth highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlolohi[i][idx[k]] has the fourth highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhihilo has the fourth lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhihilo[i][dim]
 *             has the fourth lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhihilo[i][idx[k]] has the fourth lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohilo has the third lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimlohilo[i][dim]
 *             has the third lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohilo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilolo has the second lowest doubles of the imag parts of the
 *             evaluated and differentiated monomials, outputimhilolo[i][dim]
 *             has the second lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilolo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlololo has the lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlololo[i][dim]
 *             has the lowest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimlololo[i][idx[k]] has the lowest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_dbl8_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **inputhihihi, double **inputlohihi, 
   double **inputhilohi, double **inputlolohi, 
   double **inputhihilo, double **inputlohilo, 
   double **inputhilolo, double **inputlololo, 
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cffhihihi cffhihihi[i][j] is the highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflohihi cfflohihi[i][j] is the second highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cffhilohi cffhilohi[i][j] is the third highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflolohi cfflolohi[i][j] is the fourth highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cffhihilo cffhihilo[i][j] is the fourth lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflohilo cfflohilo[i][j] is the third lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cffhilolo cffhilolo[i][j] is the second lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflololo cfflololo[i][j] is the lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   inputhihihi are the highest double coefficients of the input;
 *   inputlohihi are the second highest double coefficients of the input;
 *   inputhilohi are the third highest double coefficients of the input;
 *   inputlolohi are the fourth highest double coefficients of the input;
 *   inputhihilo are the fourth lowest double coefficients of the input;
 *   inputlohilo are the third lowest double coefficients of the input;
 *   inputhilolo are the second lowest double coefficients of the input;
 *   inputlololo are the lowest double coefficients of the input;
 *   outputhihihi has space for the output;
 *   outputlohihi has space for the output;
 *   outputhilohi has space for the output;
 *   outputlolohi has space for the output;
 *   outputhihilo has space for the output;
 *   outputlohilo has space for the output;
 *   outputhilolo has space for the output;
 *   outputlololo has space for the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputhihihi is used as work space for one column;
 *   outputlohihi is used as work space for one column;
 *   outputhilohi is used as work space for one column;
 *   outputlolohi is used as work space for one column;
 *   outputhihilo is used as work space for one column;
 *   outputlohilo is used as work space for one column;
 *   outputhilolo is used as work space for one column;
 *   outputlololo is used as work space for one column;
 *   funvalhihihi has the highest doubles of the evaluated series;
 *   funvallohihi has the second highest doubles of the evaluated series;
 *   funvalhilohi has the third highest doubles of the evaluated series;
 *   funvallolohi has the fourth highest doubles of the evaluated series;
 *   funvalhihilo has the fourth lowest doubles of the evaluated series;
 *   funvallohilo has the third lowest doubles of the evaluated series;
 *   funvalhilolo has the second lowest doubles of the evaluated series;
 *   funvallololo has the lowest doubles of the evaluated series;
 *   jacvalhihihi has the highest doubles of all derivatives;
 *   jacvallohihi has the second highest doubles of all derivatives;
 *   jacvalhilohi has the third highest doubles of all derivatives;
 *   jacvallolohi has the fourth highest doubles of all derivatives;
 *   jacvalhihilo has the fourth lowest doubles of all derivatives;
 *   jacvallohilo has the third lowest doubles of all derivatives;
 *   jacvalhilolo has the second lowest doubles of all derivatives;
 *   jacvallololo has the lowest doubles of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx8_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi, 
   double ***cffimhilohi, double ***cffimlolohi, 
   double ***cffimhihilo, double ***cffimlohilo, 
   double ***cffimhilolo, double ***cffimlololo, 
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi, 
   double ***outputimhihilo, double ***outputimlohilo, 
   double ***outputimhilolo, double ***outputimlololo, 
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cffrehihihi cffrehihihi[i][j] are the highest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelohihi cffrelohihi[i][j] are the 2nd highest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrehilohi cffrehilohi[i][j] are the 3rd highest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelolohi cffrelolohi[i][j] are the 4thd highest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrehihilo cffrehihilo[i][j] are the 4th lowest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelohilo cffrelohilo[i][j] are the 3rd lowest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrehilolo cffrehilolo[i][j] are the 2nd lowest doubles of real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelololo cffrelolo[i][j] are the lowest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhihihi cffimhihihi[i][j] are the highest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlohihi cffimlohihi[i][j] are the 2nd highest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhilohi cffimhilohi[i][j] are the 3rd highest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlolohi cffimlolohi[i][j] are the 4th highest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhihilo cffimhihilo[i][j] are the 4th lowest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlohilo cffimlohilo[i][j] are the 3rd lowest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhilolo cffimhilolo[i][j] are the 2nd lowest doubles of imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlololo cffimlololo[i][j] are the lowest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   inputrehihihi are the highest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelohihi are the 2nd highest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrehilohi are the 3rd highest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelolohi are the 4th highest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrehihilo are the 4th lowest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelohilo are the 3rd lowest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrehilolo are the 2nd lowest doubles of real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelololo are the lowest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhihihi are the highest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlohihi are the 2nd highest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhilohi are the 3rd highest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlolohi are the 4th highest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhihilo are the 4th lowest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlohilo are the 3rd lowest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhilolo are the 2nd lowest doubles of imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlololo are the lowest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   outputrehihihi has space for the highest double real parts of the output;
 *   outputrelohihi has space for the 2nd highest double real parts of output;
 *   outputrehilohi has space for the 3rd highest double real parts of output;
 *   outputrelolohi has space for the 4th highest double real parts of output;
 *   outputrehihilo has space for the 4th lowest double real parts of output;
 *   outputrelohilo has space for the 3rd lowest double real parts of output;
 *   outputrehilolo has space for the 2nd lowest double real parts of output;
 *   outputrelololo has space for the lowest double real parts of the output;
 *   outputimhihihi has space for the highest double imag parts of the output;
 *   outputimlohihi has space for the 2nd highest double imag parts of output;
 *   outputimhilohi has space for the 3rd highest double imag parts of output;
 *   outputimlolohi has space for the 4th highest double imag parts of output;
 *   outputimhihilo has space for the 4th lowest double imag parts of output;
 *   outputimlohilo has space for the 3rd lowest double imag parts of output;
 *   outputimhilolo has space for the 2nd lowest double imag parts of output;
 *   outputimlololo has space for the lowest double imag parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihihi are the highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelohihi are the 2nd highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrehilohi are the 3rd highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelolohi are the 4th highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrehihilo are the 4th lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelohilo are the 3rd lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrehilolo are the 2nd lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelololo are lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimhihihi are the highest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimlohihi are the 2nd highest double imaginary parts of evaluated
 *             and differentiated monomials used as work space for one column;
 *   outputimlohihi are the 3rd highest double imaginary parts of evaluated
 *             and differentiated monomials used as work space for one column;
 *   outputimlolohi are the 4th highest double imaginary parts of evaluated
 *             and differentiated monomials used as work space for one column;
 *   outputimhihilo are the 4th lowest double imaginary parts of evaluated
 *             differentiated monomials used as work space for one column;
 *   outputimlohilo are the 3rd lowest double imaginary parts of evaluated
 *             and differentiated monomials used as work space for one column;
 *   outputimhilolo are the 2nd lowest double imaginary parts of evaluated
 *             and differentiated monomials used as work space for one column;
 *   outputimlololo are the lowest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   funvalrehihihi are the highest double real parts of evaluated series;
 *   funvalrelohihi are the 2nd highest double real parts of evaluated series;
 *   funvalrehilohi are the 3rd highest double real parts of evaluated series;
 *   funvalrelolohi are the 4th highest double real parts of evaluated series;
 *   funvalrehihilo are the 4th lowest double real parts of evaluated series;
 *   funvalrelohilo are the 3rd lowest double real parts of evaluated series;
 *   funvalrehilolo are the 2nd lowest double real parts of evaluated series;
 *   funvalrelololo are the lowest double real parts of evaluated series;
 *   funvalimhihihi are the highest double imag parts of evaluated series;
 *   funvalimlohihi are the 2nd highest double imag parts of evaluated series;
 *   funvalimhilohi are the 3rd highest double imag parts of evaluated series;
 *   funvalimlolohi are the 4th highest double imag parts of evaluated series;
 *   funvalimhihilo are the 4th lowest double imag parts of evaluated series;
 *   funvalimlohilo are the 3rd lowest double imag parts of evaluated series;
 *   funvalimhilolo are the 2nd lowest double imag parts of evaluated series;
 *   funvalimlololo are the lowest double imag parts of evaluated series;
 *   jacvalrehihihi are the highest double real parts of all derivatives;
 *   jacvalrelohihi are the 2nd highest double real parts of all derivatives;
 *   jacvalrehilohi are the 3rd highest double real parts of all derivatives;
 *   jacvalrelolhi are the 4th highest double real parts of all derivatives;
 *   jacvalrehihilo are the 4th lowest double real parts of all derivatives;
 *   jacvalrelohilo are the 3rd lowest double real parts of all derivatives;
 *   jacvalrehilolo are the 2nd lowest double real parts of all derivatives;
 *   jacvalrelololo are the lowest double real parts of all derivatives;
 *   jacvalimhihihi are the highest double imag parts of all derivatives;
 *   jacvalimlohihi are the 2nd highest double imag parts of all derivatives;
 *   jacvalimhilohi are the 3rd highest double imag parts of all derivatives;
 *   jacvalimlolohi are the 4th highest double imag parts of all derivatives;
 *   jacvalimhihilo are the 4th lowest double imag parts of all derivatives;
 *   jacvalimlohilo are the 3rd lowest double imag parts of all derivatives;
 *   jacvalimhilolo are the 2nd lowest double imag parts of all derivatives;
 *   jacvalimlololo are the lowest double imag parts of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

#endif
