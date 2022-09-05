// The file dbl8_systems_kernels.h specifies functions to define
// memory transfers and kernel launches to evaluate and differentiate
// monomials with common factors in octo double precision.

#ifndef __dbl8_systems_kernels_h__
#define __dbl8_systems_kernels_h__

#include "convolution_jobs.h"

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
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose );
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
 *   verbose  if true, writes extra information.
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
   double *cnvlapms, double *elapsedms, double *walltimesec, bool verbose );
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
 *   verbose  if true, then extra output about the setup is written.
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
   double ***outputhilolo, double ***outputlololo, int vrblvl );
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
 *   vrblvl    is the verbose level.
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
 *            outputlololo[dim] has the lowest double value. */

#endif
