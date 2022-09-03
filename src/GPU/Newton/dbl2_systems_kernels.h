// The file dbl2_systems_kernels.h specifies functions to define
// memory transfers ad kernel launches to evaluate and differentiate
// monomials with common factors in double double precision.

#ifndef __dbl2_systems_kernels_h__
#define __dbl2_systems_kernels_h__

#include "convolution_jobs.h"

void dbl2_evaldiffdata_to_output
 ( double *datahi, double *datalo, double ***outputhi, double ***outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *
 * ON ENTRY :
 *   datahi   high doubles of coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products;
 *   datalo   low doubles of coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products;
 *   outputhi has space for the high double values and all derivatives;
 *   outputlo has space for the low double values and all derivatives;
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
 *   outputhi has in outputhi[i][dim] the highest double value
 *            of the i-th monomial, and
 *            outputhi[i][k] has the highest double value of
 *            the k-th derivative of the i-th monomial;
 *   outputlo has in outputlo[i][dim] the lowest double value
 *            of the i-th monomial, and
 *            outputlo[i][k] has the lowest double value of
 *            the k-th derivative of the i-th monomial. */

void GPU_dbl2_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhi, double **cfflo, double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo, ConvolutionJobs cnvjobs,
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
 *   cffhi    cffhi[k] has deg+1 doubles for the high doubles of 
 *            the coefficient series of monomial k;
 *   cfflo    cfflo[k] has deg+1 doubles for the low doubles of 
 *            the coefficient series of monomial k;
 *   inputhi  has the high doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputlo  has the low doubles of the coefficients of the series
 *            for all variables in the polynomial;
 *   outputhi has space allocated for the high double values and derivatives;
 *   outputlo has space allocated for the low double values and derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   verbose  if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputhi has the high doubles of the output,
 *            outputhi[k], for k from 0 to dim-1, has the high double of
 *            the derivative with respect to the variable k;
 *            outputhi[dim] has the high double value of the polynomial;
 *   outputlo has the low doubles of the output,
 *            outputlo[k], for k from 0 to dim-1, has the low double of
 *            the derivative with respect to the variable k;
 *            outputlo[dim] has the low double value of the polynomial;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_dbl2_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo, double ***outputhi, double ***outputlo,
   int vrblvl );
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
 *   cffhi    high doubles of the coefficients of the monomials;
 *   cfflo    low doubles of the coefficients of the monomials;
 *   acchi    space to accumulate one power series of degree deg;
 *   acclo    space to accumulate one power series of degree deg;
 *   inputhi  high doubles of the coefficients of the series of degree deg,
 *            for dim variables;
 *   inputlo  low doubles of the coefficients of the series of degree deg,
 *            for dim variables;
 *   outputhi has space for the high doubles of the output;
 *   outputlo has space for the low doubles of the output;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cffhi    has the high doubles of the evaluated common factors;
 *   cfflo    has the low doubles of the evaluated common factors;
 *   outputhi has the high doubles of the output,
 *            outputhi[k], for k from 0 to dim-1, has the high double of
 *            the derivative with respect to the variable k;
 *            outputhi[dim] has the high double value of the polynomial;
 *   outputlo has the low doubles of the output,
 *            outputlo[k], for k from 0 to dim-1, has the low double of
 *            the derivative with respect to the variable k;
 *            outputlo[dim] has the low double value of the polynomial. */

#endif
