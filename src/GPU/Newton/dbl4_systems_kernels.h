// The file dbl4_systems_kernels.h specifies functions to define
// memory transfers ad kernel launches to evaluate and differentiate
// monomials with common factors in quad double precision.

#ifndef __dbl4_systems_kernels_h__
#define __dbl4_systems_kernels_h__

#include "convolution_jobs.h"

void dbl4_evaldiffdata_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose );
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
 *   verbose  if true, writes extra information.
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

void GPU_dbl4_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, ConvolutionJobs cnvjobs,
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
 *   verbose  if true, then extra output about the setup is written.
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

void GPU_dbl4_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, int vrblvl );
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
 *   vrblvl    is the verbose level.
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
 *            outputlolo[dim] has the lowest double value of the polynomial. */

#endif
