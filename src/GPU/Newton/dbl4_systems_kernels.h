// The file dbl4_systems_kernels.h specifies functions to define
// memory transfers ad kernel launches to evaluate and differentiate
// monomials with common factors in quad double precision.

#ifndef __dbl4_systems_kernels_h__
#define __dbl4_systems_kernels_h__

#include "convolution_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"

void write_dbl4_cnvflops
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

void write_vectorized4_cnvincflops
 ( int dim, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the kernel and wall time flops for the convolution
 *   and the increment jobs in the complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim      number of series;
 *   deg      order of the series, truncation degree;
 *   cnvjobs  defines the convolution jobs;
 *   incjobs  defines the increment jobs;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void dbl4_evaldiffdata_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *
 * ON ENTRY :
 *   datahihi has the highest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalohi has second highest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   datahilo has second lowest doubles of monomials and input, 
 *            computed forward, backward, and cross products;
 *   datalolo has the lowest doubles of monomials and input,
 *            computed forward, backward, and cross products;
 *   outputhihi has space for the highest double values and derivatives;
 *   outputlohi has space for the second highest double values and derivatives;
 *   outputhilo has space for the second lowest double values and derivatives;
 *   outputlolo has space for the lowest double values and derivatives;
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
 *   outputhihi has in outputhihi[i][dim] the highest double value
 *            of the i-th monomial, and
 *            outputhi[i][k] has the highest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputlohi has in outputlohi[i][dim] the second highest double value
 *            of the i-th monomial, and
 *            outputlo[i][k] has the second highest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputhilo has in outputhilo[i][dim] the second lowest double value
 *            of the i-th monomial, and
 *            outputhilo[i][k] has the second lowest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputlolo has in outputlolo[i][dim] the lowest double value
 *            of the i-th monomial, and
 *            outputlolo[i][k] has the lowest double value of
 *            the k-th derivative of the i-th monomial. */

void cmplx4_evaldiffdata_to_output
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *
 * ON ENTRY :
 *   datarehihi are the highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelohi are the second highest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarehilo are the second lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   datarelolo are the lowest doubles of the real parts of the input, 
 *            computed forward, backward, and cross products;
 *   dataimhihi are the highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlohi are the second highest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimhilo are the second lowest doubles of the imaginary parts of the
 *            input, computed forward, backward, and cross products;
 *   dataimlolo are the lowest doubles of the imaginary parts of the 
 *            input, computed forward, backward, and cross products;
 *   outputrehihi has space for the value and all derivatives;
 *   outputrelohi has space for the value and all derivatives;
 *   outputrehilo has space for the value and all derivatives;
 *   outputrelolo has space for the value and all derivatives;
 *   outputimhihi has space for the value and all derivatives;
 *   outputimlohi has space for the value and all derivatives;
 *   outputimhilo has space for the value and all derivatives;
 *   outputimlolo has space for the value and all derivatives;
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
 *   outputrehihi, in outputrehihi[i][dim] are the highest doubles of the real
 *            parts of the value of the i-th monomial, and outputrehihi[i][k]
 *            has the highest doubles of the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputrelohi, in outputrelohi[i][dim] are the second highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelohi[i][k] has the second highest doubles of the
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilo, in outputrehilo[i][dim] are the second lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehilo[i][k] has the lowest doubles of the real parts of
 *            the value of the k-th derivative of the i-th monomial;
 *   outputrelolo, in outputrelolo[i][dim] are the lowest doubles of the real
 *            parts of the value of the i-th monomial, and outputrelolo[i][k]
 *            has the lowest doubles of the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputimhihi, in outputrehihi[i][dim] are the highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhihi[i][k] has the highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohi, in outputrelohi[i][dim] are the second highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohi[i][k] has the second highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilo, in outputrehilo[i][dim] are the second lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohi[i][k] has the second lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlolo, in outputrelolo[i][dim] are the lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolo[i][k] has the lowest doubles of the imaginary parts
 *            of the value of the k-th derivative of the i-th monomial. */

void cmplx4vectorized_evaldiffdata_to_output
 ( double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   int dim, int nbr, int deg, int *nvr, int totcffoffset,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device,
 *   using vectorized complex arithmetic, to the output.
 *
 * ON ENTRY :
 *   datarihihi are the highest doubles of the real/imag parts
 *            of the input, computed forward, backward, and cross products;
 *   datarilohi are the second highest doubles of the real/imag part 
 *            of the input, computed forward, backward, and cross products;
 *   datarihilo are the second lowest doubles of the real/imag parts
 *            of the input, computed forward, backward, and cross products;
 *   datarilolo are the lowest doubles of the real/imag parts 
 *            of the input, computed forward, backward, and cross products;
 *   outputrehihi has space for the value and all derivatives;
 *   outputrelohi has space for the value and all derivatives;
 *   outputrehilo has space for the value and all derivatives;
 *   outputrelolo has space for the value and all derivatives;
 *   outputimhihi has space for the value and all derivatives;
 *   outputimlohi has space for the value and all derivatives;
 *   outputimhilo has space for the value and all derivatives;
 *   outputimlolo has space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   totcffoffset is the offset for the imaginary parts of the data;
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
 *   outputrehihi, in outputrehihi[i][dim] are the highest doubles of the real
 *            parts of the value of the i-th monomial, and outputrehihi[i][k]
 *            has the highest doubles of the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputrelohi, in outputrelohi[i][dim] are the second highest doubles of
 *            the real parts of the value of the i-th monomial, and
 *            outputrelohi[i][k] has the second highest doubles of the
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputrehilo, in outputrehilo[i][dim] are the second lowest doubles of
 *            the real parts of the value of the i-th monomial, and 
 *            outputrehilo[i][k] has the lowest doubles of the real parts of
 *            the value of the k-th derivative of the i-th monomial;
 *   outputrelolo, in outputrelolo[i][dim] are the lowest doubles of the real
 *            parts of the value of the i-th monomial, and outputrelolo[i][k]
 *            has the lowest doubles of the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputimhihi, in outputrehihi[i][dim] are the highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimhihi[i][k] has the highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlohi, in outputrelohi[i][dim] are the second highest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohi[i][k] has the second highest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimhilo, in outputrehilo[i][dim] are the second lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlohi[i][k] has the second lowest doubles of the imaginary
 *            parts of the value of the k-th derivative of the i-th monomial;
 *   outputimlolo, in outputrelolo[i][dim] are the lowest doubles of
 *            the imaginary parts of the value of the i-th monomial, and
 *            outputimlolo[i][k] has the lowest doubles of the imaginary parts
 *            of the value of the k-th derivative of the i-th monomial. */

void GPU_dbl4_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, ConvolutionJobs cnvjobs,
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
 *   cffhihi  cffhihi[k] has deg+1 doubles for the highest doubles
 *            of the coefficient series of monomial k;
 *   cfflohi  cfflohi[k] has deg+1 doubles for the second highest doubles
 *            of the coefficient series of monomial k;
 *   cffhilo  cffhilo[k] has deg+1 doubles for the second lowest doubles
 *            of the coefficient series of monomial k;
 *   cfflolo  cfflolo[k] has deg+1 doubles for the lowest doubles
 *            of the coefficient series of monomial k;
 *   inputhihi has the highest doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputlohi has the second highest doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputhilo has the second lowest doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputlolo has the lowest doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   outputhihi has space allocated for the highest double output;
 *   outputlohi has space allocated for the second highest double output;
 *   outputhilo has space allocated for the second lowest double output;
 *   outputlolo has space allocated for the lowest double output;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputhihi has the highest doubles of the output,
 *            outputhihi[k], for k from 0 to dim-1, has the highest double
 *            of the derivative with respect to the variable k;
 *            outputhihi[dim] has the highest double value of the polynomial;
 *   outputlohi has the second highest doubles of the output,
 *            outputlohi[k], for k from 0 to dim-1, has the second highest
 *            double of the derivative with respect to the variable k;
 *            outputlohi[dim] has the seconc highest double value;
 *   outputhilo has the second lowest doubles of the output,
 *            outputhilo[k], for k from 0 to dim-1, has the second lowest
 *            double of the derivative with respect to the variable k;
 *            outputhilo[dim] has the second lowest double value;
 *   outputlolo has the lowest doubles of the output,
 *            outputlolo[k], for k from 0 to dim-1, has the lowest double
 *            of the derivative with respect to the variable k;
 *            outputlolo[dim] has the lowest double value of the polynomial;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx4_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo, ConvolutionJobs cnvjobs,
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
 *   cffrehihi, cffrehihi[k] has deg+1 highest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffrelohi, cffrelohi[k] has deg+1 second highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilo, cffrehilo[k] has deg+1 second lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelolo, cffrelolo[k] has deg+1 lowest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffimhihi, cffimhihi[k] has deg+1 highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohi, cffimlohi[k] has deg+1 second highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilo, cffimihilo[k] has deg+1 second lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlolo, cffimlolo[k] has deg+1 lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   inputrehihi are the highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohi are the second highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilo are the second lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelolo are the lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihi are the highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohi are the second highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilo are the second lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlolo are the lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   outputrehihi has space for the highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohi has space for the second highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilo has space for the second lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelolo has space for the lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputimhihi has space for the highest doubles of the imaginary parts
 *            of value and all derivatives;
 *   outputimlohi has space for the second highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilo has space for the second lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlolo has space for the lowest doubles of the imaginary parts
 *            of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihi has the highest doubles of the real parts of derivatives
 *            and the value, outputrehihi[k], for k from 0 to dim-1,
 *            has the highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihi[dim] has the
 *            highest doubles of the real part of the value;
 *   outputrelohi has the 2nd highest doubles of the real parts of derivatives
 *            and the value, outputrelohi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohi[dim] has the
 *            second highest doubles of the real part of the value;
 *   outputrehilo has the 2nd lowest doubles of the real parts of derivatives
 *            and the value, outputrehilo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilo[dim] has the
 *            second lowest doubles of the real part of the value;
 *   outputrelolo has the lowest doubles of the real parts of derivatives
 *            and the value, outputrelolo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelolo[dim] has the
 *            lowest doubles of the real part of the value;
 *   outputimhihi has the highest doubles of the imaginary parts of derivatives
 *            and the value, outputimhihi[k], for k from 0 to dim-1, has the
 *            highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihi[dim] has the
 *            highest doubles of the imaginary part of the value;
 *   outputimlohi has the 2nd highest doubles of the imag parts of derivatives
 *            and the value, outputimlohi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohi[dim] has the
 *            second highest doubles of the imaginary part of the value;
 *   outputimhilo has the 2nd lowest doubles of the imag parts of derivatives
 *            and the value, outputimhilo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilo[dim] has the
 *            second lowest doubles of the imaginary part of the value;
 *   outputimlolo has the lowest doubles of the imaginary parts of derivatives
 *            and the value, outputimlolo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlolo[dim] has the
 *            lowest doubles of the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx4vectorized_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   with vectorized complex arithmetic.
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
 *   cffrehihi, cffrehihi[k] has deg+1 highest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffrelohi, cffrelohi[k] has deg+1 second highest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrehilo, cffrehilo[k] has deg+1 second lowest doubles of the real
 *            parts for the coefficient series of monomial k;
 *   cffrelolo, cffrelolo[k] has deg+1 lowest doubles of the real parts
 *            for the coefficient series of monomial k;
 *   cffimhihi, cffimhihi[k] has deg+1 highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlohi, cffimlohi[k] has deg+1 second highest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimhilo, cffimihilo[k] has deg+1 second lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   cffimlolo, cffimlolo[k] has deg+1 lowest doubles of the
 *            imaginary parts for the coefficient series of monomial k;
 *   inputrehihi are the highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelohi are the second highest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrehilo are the second lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputrelolo are the lowest doubles of the real parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhihi are the highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlohi are the second highest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimhilo are the second lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   inputimlolo are the lowest doubles of the imaginary parts of the
 *            coefficients of the series for all variables in the polynomial;
 *   outputrehihi has space for the highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelohi has space for the second highest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrehilo has space for the second lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputrelolo has space for the lowest doubles of the real parts
 *            of the value and all derivatives;
 *   outputimhihi has space for the highest doubles of the imaginary parts
 *            of value and all derivatives;
 *   outputimlohi has space for the second highest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimhilo has space for the second lowest doubles of the imaginary
 *            parts of value and all derivatives;
 *   outputimlolo has space for the lowest doubles of the imaginary parts
 *            of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   incjobs  increment jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihi has the highest doubles of the real parts of derivatives
 *            and the value, outputrehihi[k], for k from 0 to dim-1,
 *            has the highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehihi[dim] has the
 *            highest doubles of the real part of the value;
 *   outputrelohi has the 2nd highest doubles of the real parts of derivatives
 *            and the value, outputrelohi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelohi[dim] has the
 *            second highest doubles of the real part of the value;
 *   outputrehilo has the 2nd lowest doubles of the real parts of derivatives
 *            and the value, outputrehilo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrehilo[dim] has the
 *            second lowest doubles of the real part of the value;
 *   outputrelolo has the lowest doubles of the real parts of derivatives
 *            and the value, outputrelolo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the real part of the derivative
 *            with respect to the variable k; outputrelolo[dim] has the
 *            lowest doubles of the real part of the value;
 *   outputimhihi has the highest doubles of the imaginary parts of derivatives
 *            and the value, outputimhihi[k], for k from 0 to dim-1, has the
 *            highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhihi[dim] has the
 *            highest doubles of the imaginary part of the value;
 *   outputimlohi has the 2nd highest doubles of the imag parts of derivatives
 *            and the value, outputimlohi[k], for k from 0 to dim-1, has the
 *            second highest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlohi[dim] has the
 *            second highest doubles of the imaginary part of the value;
 *   outputimhilo has the 2nd lowest doubles of the imag parts of derivatives
 *            and the value, outputimhilo[k], for k from 0 to dim-1, has the
 *            second lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimhilo[dim] has the
 *            second lowest doubles of the imaginary part of the value;
 *   outputimlolo has the lowest doubles of the imaginary parts of derivatives
 *            and the value, outputimlolo[k], for k from 0 to dim-1, has the
 *            lowest doubles of the imaginary part of the derivative
 *            with respect to the variable k; outputimlolo[dim] has the
 *            lowest doubles of the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_dbl4_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
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
 *   cffhihi  highest doubles of the coefficients of the monomials;
 *   cfflohi  second highest doubles of the coefficients of the monomials;
 *   cffhilo  second lowest doubles of the coefficients of the monomials;
 *   cfflolo  lowest doubles of the coefficients of the monomials;
 *   acchihi  space to accumulate one power series of degree deg;
 *   acclolo  space to accumulate one power series of degree deg;
 *   acchihi  space to accumulate one power series of degree deg;
 *   acclolo  space to accumulate one power series of degree deg;
 *   inputhihi has the highest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputlohi has the second highest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputhilo has the second lowest doubles of the coefficients of the input,
 *            for dim variables;
 *   inputlolo has the lowest doubles of the coefficients of the input,
 *            for dim variables;
 *   outputhihi has space for the highest doubles of the output;
 *   outputlohi has space for the second highest doubles of the output;
 *   outputhilo has space for the second lowest doubles of the output;
 *   outputlolo has space for the lowest doubles of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffhihi  highest doubles of the evaluated common factors;
 *   cfflohi  second highest doubles of the evaluated common factors;
 *   cffhilo  second lowest doubles of the evaluated common factors;
 *   cfflolo  lowest doubles of the evaluated common factors;
 *   outputhihi has the highest doubles of the output,
 *            outputhihi[k], for k from 0 to dim-1, has the highest double
 *            of the derivative with respect to the variable k;
 *            outputhihi[dim] has the highest double value of the polynomial;
 *   outputlohi has the second highest doubles of the output,
 *            outputlohi[k], for k from 0 to dim-1, has the second highest
 *            double of the derivative with respect to the variable k;
 *            outputlohi[dim] has the seconc highest double value;
 *   outputhilo has the second lowest doubles of the output,
 *            outputhilo[k], for k from 0 to dim-1, has the second lowest
 *            double of the derivative with respect to the variable k;
 *            outputhilo[dim] has the second lowest double value;
 *   outputlolo has the lowest doubles of the output,
 *            outputlolo[k], for k from 0 to dim-1, has the lowest double
 *            of the derivative with respect to the variable k;
 *            outputlolo[dim] has the lowest double value of the polynomial;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx4_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi,
   double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi,
   double *accimhilo, double *accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi, 
   double **inputimhilo, double **inputimlolo, 
   double ***outputrehihi, double ***outputrelohi, 
   double ***outputrehilo, double ***outputrelolo, 
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
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
 *   cffrehihi are the highest doubles of the real parts of the coefficients;
 *   cffrelohi are the second highest doubles of the real parts
 *             of the coefficients;
 *   cffrehilo are the second lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelolo are the lowest doubles of the real parts of the coefficients
 *   cffimhihi are the highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohi are the second highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilo are the second lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlolo are the lowest doubles of the imaginary parts
 *             of the coefficients;
 *   accrehihi has space to accumulate one series of degree deg;
 *   accrelohi has space to accumulate one series of degree deg;
 *   accrehilo has space to accumulate one series of degree deg;
 *   accrelolo has space to accumulate one series of degree deg;
 *   accimhihi has space to accumulate one  series of degree deg;
 *   accimlohi has space to accumulate one  series of degree deg;
 *   accimhilo has space to accumulate one  series of degree deg;
 *   accimlolo has space to accumulate one  series of degree deg;
 *   inputrehihi has the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi has the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo has the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelolo has the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihi has the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi has the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo has the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolo has the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihi has space for the highest doubles of the real parts
 *             of the output;
 *   outputrelohi has space for the second highest doubles of the real parts
 *             of the output;
 *   outputrehilo has space for the second lowest doubles of the real parts
 *             of the output;
 *   outputrelolo has space for the lowest doubles of the real parts
 *             of the output;
 *   outputimhihi has space for the highest doubles of the imaginary
 *             parts of the output;
 *   outputimlohi has space for the second highest doubles of the imaginary
 *             parts of the output;
 *   outputimhilo has space for the second lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlolo has space for the lowest doubles of the imaginary
 *             parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffrehihi has the highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohi has the second highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilo has the second lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelolo has the lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffimhihi has the highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohi has the second highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilo has the second lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlolo has the lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   outputrehihi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrehihi[i][dim] has the
 *             highest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihi[i][idx[k]] has the highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelohi[i][dim] has the
 *             second highest doubles of the real parts of the value of the
 *             input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohi[i][idx[k]] has the second highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilo has the second lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilo[i][dim]
 *             has the second lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilo[i][idx[k]] has the second lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelolo has the lowest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelolo[i][dim] has the
 *             lowest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelolo[i][idx[k]] has the lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputimhihi has the highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhihi[i][dim]
 *             has the highest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimhi[i][idx[k]] has the highest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohi has the second highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlohi[i][dim]
 *             has the second highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohi[i][idx[k]] has the second highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilo has the second lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhilo[i][dim]
 *             has the second lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlolo has the lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlolo[i][dim]
 *             has the lowest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimlolo[i][idx[k]] has the lowest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx4vectorized_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi,
   double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi,
   double *accimhilo, double *accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi, 
   double **inputimhilo, double **inputimlolo, 
   double ***outputrehihi, double ***outputrelohi, 
   double ***outputrehilo, double ***outputrelolo, 
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series using complex vectorized arithmetic.
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
 *   cffrehihi are the highest doubles of the real parts of the coefficients;
 *   cffrelohi are the second highest doubles of the real parts
 *             of the coefficients;
 *   cffrehilo are the second lowest doubles of the real parts
 *             of the coefficients;
 *   cffrelolo are the lowest doubles of the real parts of the coefficients
 *   cffimhihi are the highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlohi are the second highest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimhilo are the second lowest doubles of the imaginary parts
 *             of the coefficients;
 *   cffimlolo are the lowest doubles of the imaginary parts
 *             of the coefficients;
 *   accrehihi has space to accumulate one series of degree deg;
 *   accrelohi has space to accumulate one series of degree deg;
 *   accrehilo has space to accumulate one series of degree deg;
 *   accrelolo has space to accumulate one series of degree deg;
 *   accimhihi has space to accumulate one  series of degree deg;
 *   accimlohi has space to accumulate one  series of degree deg;
 *   accimhilo has space to accumulate one  series of degree deg;
 *   accimlolo has space to accumulate one  series of degree deg;
 *   inputrehihi has the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi has the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo has the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelolo has the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihi has the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi has the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo has the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolo has the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihi has space for the highest doubles of the real parts
 *             of the output;
 *   outputrelohi has space for the second highest doubles of the real parts
 *             of the output;
 *   outputrehilo has space for the second lowest doubles of the real parts
 *             of the output;
 *   outputrelolo has space for the lowest doubles of the real parts
 *             of the output;
 *   outputimhihi has space for the highest doubles of the imaginary
 *             parts of the output;
 *   outputimlohi has space for the second highest doubles of the imaginary
 *             parts of the output;
 *   outputimhilo has space for the second lowest doubles of the imaginary
 *             parts of the output;
 *   outputimlolo has space for the lowest doubles of the imaginary
 *             parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cffrehihi has the highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelohi has the second highest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrehilo has the second lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffrelolo has the lowest doubles of the real parts
 *             of the evaluated common factors;
 *   cffimhihi has the highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlohi has the second highest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimhilo has the second lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   cffimlolo has the lowest doubles of the imaginary parts
 *             of the evaluated common factors;
 *   outputrehihi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrehihi[i][dim] has the
 *             highest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehihi[i][idx[k]] has the highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelohi has the highest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelohi[i][dim] has the
 *             second highest doubles of the real parts of the value of the
 *             input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelohi[i][idx[k]] has the second highest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrehilo has the second lowest doubles of the real parts of the
 *             evaluated and differentiated monomials, outputrehilo[i][dim]
 *             has the second lowest doubles of the real parts of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrehilo[i][idx[k]] has the second lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputrelolo has the lowest doubles of the real parts of the evaluated
 *             and differentiated monomials, outputrelolo[i][dim] has the
 *             lowest doubles of the real parts of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputrelolo[i][idx[k]] has the lowest doubles
 *             of the real parts of the derivative w.r.t. idx[k];
 *   outputimhihi has the highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhihi[i][dim]
 *             has the highest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimhi[i][idx[k]] has the highest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlohi has the second highest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlohi[i][dim]
 *             has the second highest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimlohi[i][idx[k]] has the second highest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimhilo has the second lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimhilo[i][dim]
 *             has the second lowest doubles of the imaginary part of the
 *             value of the input at the i-th monomial, and for k in range
 *             0..nvr[i]-1: outputimhilo[i][idx[k]] has the second lowest
 *             doubles of the imaginary parts of te derivative w.r.t. idx[k];
 *   outputimlolo has the lowest doubles of the imaginary parts of the
 *             evaluated and differentiated monomials, outputimlolo[i][dim]
 *             has the lowest doubles of the imaginary part of the value of
 *             the input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputimlolo[i][idx[k]] has the lowest doubles
 *             of the imaginary parts of te derivative w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_dbl4_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **inputhihi, double **inputlohi, 
   double **inputhilo, double **inputlolo, 
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo,
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
 *   cffhihi   cffhihi[i][j] is the highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflohi   cfflohi[i][j] is the second highest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cffhilo   cffhilo[i][j] is the second lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   cfflolo   cfflolo[i][j] is the lowest double of the coefficient
 *             of the j-th monomial in the i-th column;
 *   inputhihi are the highest double coefficients of the input;
 *   inputlohi are the second highest double coefficients of the input;
 *   inputhilo are the second lowest double coefficients of the input;
 *   inputlolo are the lowest double coefficients of the input;
 *   outputhihi has space for the output;
 *   outputlohi has space for the output;
 *   outputhilo has space for the output;
 *   outputlolo has space for the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputhihi is used as work space for one column;
 *   outputlohi is used as work space for one column;
 *   outputhilo is used as work space for one column;
 *   outputlolo is used as work space for one column;
 *   funvalhihi has the highest doubles of the evaluated series;
 *   funvallohi has the second highest doubles of the evaluated series;
 *   funvalhilo has the second lowest doubles of the evaluated series;
 *   funvallolo has the lowest doubles of the evaluated series;
 *   jacvalhihi has the highest doubles of all derivatives;
 *   jacvallohi has the second highest doubles of all derivatives;
 *   jacvalhilo has the second lowest doubles of all derivatives;
 *   jacvallolo has the lowest doubles of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx4_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi, 
   double ***cffimhilo, double ***cffimlolo, 
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi, 
   double ***outputimhilo, double ***outputimlolo, 
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo,
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
 *   cffrehihi cffrehihi[i][j] are the highest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelohi cffrelohi[i][j] are the 2nd highest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrehilo cffrehilo[i][j] are the 2nd lowest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelolo cffrelolo[i][j] are the lowest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhihi cffimhihi[i][j] are the highest doubles of the imaginary parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlohi cffimlohi[i][j] are the 2nd highest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhilo cffimhilo[i][j] are the 2nd lowest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlolo cffimlolo[i][j] are the lowest doubles of the imaginary parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   inputrehihi are the highest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelohi are the 2nd highest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrehilo are the 2nd lowest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelolo are the lowest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhihi are the highest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlohi are the 2nd highest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhilo are the 2nd lowest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlolo are the lowest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   outputrehihi has space for the highest double real parts of the output;
 *   outputrelohi has space for the 2nd highest double real parts of output;
 *   outputrehilo has space for the 2nd lowest double real parts of output;
 *   outputrelolo has space for the lowest double real parts of the output;
 *   outputimhihi has space for the highest double imag parts of the output;
 *   outputimlohi has space for the 2nd highest double imag parts of output;
 *   outputimhilo has space for the 2nd lowest double imag parts of output;
 *   outputimlolo has space for the lowest double imag parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihi are the highiest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelohi are the 2nd highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrehilo are the 2nd lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelolo are lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimhihi are the highest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimlohi are the 2nd highest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimhilo are the 2nd lowest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimlolo are the lowest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   funvalrehihi are the highest double real parts of evaluated series;
 *   funvalrelohi are the 2nd highest double real parts of evaluated series;
 *   funvalrehilo are the 2nd lowest double real parts of evaluated series;
 *   funvalrelolo are the lowest double real parts of evaluated series;
 *   funvalimhihi are the highest double imag parts of evaluated series;
 *   funvalimlohi are the 2nd highest double imag parts of evaluated series;
 *   funvalimhilo are the 2nd lowest double imag parts of evaluated series;
 *   funvalimlolo are the lowest double imag parts of evaluated series;
 *   jacvalrehihi are the highest double real parts of all derivatives;
 *   jacvalrelohi are the 2nd highest double real parts of all derivatives;
 *   jacvalrehilo are the 2nd lowest double real parts of all derivatives;
 *   jacvalrelolo are the lowest double real parts of all derivatives;
 *   jacvalimhihi are the highest double imag parts of all derivatives;
 *   jacvalimlohi are the 2nd highest double imag parts of all derivatives;
 *   jacvalimhilo are the 2nd lowest double imag parts of all derivatives;
 *   jacvalimlolo are the lowest double imag parts of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx4vectorized_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi, 
   double ***cffimhilo, double ***cffimlolo, 
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi, 
   double ***outputimhilo, double ***outputimlolo, 
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, with vectorized complex arithmetic.
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
 *   cffrehihi cffrehihi[i][j] are the highest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelohi cffrelohi[i][j] are the 2nd highest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrehilo cffrehilo[i][j] are the 2nd lowest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffrelolo cffrelolo[i][j] are the lowest doubles of the real parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhihi cffimhihi[i][j] are the highest doubles of the imaginary parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlohi cffimlohi[i][j] are the 2nd highest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimhilo cffimhilo[i][j] are the 2nd lowest doubles of the imag parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   cffimlolo cffimlolo[i][j] are the lowest doubles of the imaginary parts
 *             of the coefficients of the j-th monomial in the i-th column;
 *   inputrehihi are the highest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelohi are the 2nd highest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrehilo are the 2nd lowest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelolo are the lowest doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhihi are the highest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlohi are the 2nd highest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhilo are the 2nd lowest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlolo are the lowest doubles of the imag parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   outputrehihi has space for the highest double real parts of the output;
 *   outputrelohi has space for the 2nd highest double real parts of output;
 *   outputrehilo has space for the 2nd lowest double real parts of output;
 *   outputrelolo has space for the lowest double real parts of the output;
 *   outputimhihi has space for the highest double imag parts of the output;
 *   outputimlohi has space for the 2nd highest double imag parts of output;
 *   outputimhilo has space for the 2nd lowest double imag parts of output;
 *   outputimlolo has space for the lowest double imag parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputrehihi are the highiest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelohi are the 2nd highest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrehilo are the 2nd lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputrelolo are lowest double real parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimhihi are the highest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimlohi are the 2nd highest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimhilo are the 2nd lowest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   outputimlolo are the lowest double imaginary parts of evaluated and
 *             differentiated monomials used as work space for one column;
 *   funvalrehihi are the highest double real parts of evaluated series;
 *   funvalrelohi are the 2nd highest double real parts of evaluated series;
 *   funvalrehilo are the 2nd lowest double real parts of evaluated series;
 *   funvalrelolo are the lowest double real parts of evaluated series;
 *   funvalimhihi are the highest double imag parts of evaluated series;
 *   funvalimlohi are the 2nd highest double imag parts of evaluated series;
 *   funvalimhilo are the 2nd lowest double imag parts of evaluated series;
 *   funvalimlolo are the lowest double imag parts of evaluated series;
 *   jacvalrehihi are the highest double real parts of all derivatives;
 *   jacvalrelohi are the 2nd highest double real parts of all derivatives;
 *   jacvalrehilo are the 2nd lowest double real parts of all derivatives;
 *   jacvalrelolo are the lowest double real parts of all derivatives;
 *   jacvalimhihi are the highest double imag parts of all derivatives;
 *   jacvalimlohi are the 2nd highest double imag parts of all derivatives;
 *   jacvalimhilo are the 2nd lowest double imag parts of all derivatives;
 *   jacvalimlolo are the lowest double imag parts of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

#endif
