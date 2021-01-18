// The file dbl_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in double precision.

#ifndef __dbl_polynomials_kernels_h__
#define __dbl_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

__global__ void dbl_padded_convjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all convolution jobs at the same layer, on real data.
 *   The block index defines the convolution job.
 *
 * REQUIRED : 
 *   The number of blocks equals the size of  in1idx, in2idx, outidx,
 *   and dim equals the number of threads in each block.
 *
 * ON ENTRY :
 *   data      coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the convolution jobs;
 *   in2idx    indices of the second input of the convolution jobs;
 *   outidx    indices of the output of the convolution jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   data      updated forward, backward, and cross products. */

__global__ void cmplx_padded_convjobs
 ( double *datare, double *dataim,
   int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all convolution jobs at the same layer, on complex data.
 *   The block index defines the convolution job.
 *
 * REQUIRED : 
 *   The number of blocks equals the size of  in1idx, in2idx, outidx,
 *   and dim equals the number of threads in each block.
 *
 * ON ENTRY :
 *   datare    real parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   dataim    imaginary parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the convolution jobs;
 *   in2idx    indices of the second input of the convolution jobs;
 *   outidx    indices of the output of the convolution jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   datare    updated real parts of 
 *             the forward, backward, and cross products;
 *   dataim    updated imaginary parts of 
 *             the forward, backward, and cross products. */

__global__ void dbl_update_addjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all addition jobs at the same layer, on real data.
 *   The block index defines the addition job.
 *
 * REQUIRED : 
 *   The number of blocks equals the size of  in1idx, in2idx, outidx,
 *   and dim equals the number of threads in each block.
 *
 * ON ENTRY :
 *   data      coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the addition jobs;
 *   in2idx    indices of the second input of the addition jobs;
 *   outidx    indices of the output of the addition jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   data      updated forward, backward, and cross products. */

__global__ void cmplx_update_addjobs
 ( double *datare, double *dataim,
   int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all addition jobs at the same layer, on complex data.
 *   The block index defines the addition job.
 *
 * REQUIRED : 
 *   The number of blocks equals the size of  in1idx, in2idx, outidx,
 *   and dim equals the number of threads in each block.
 *
 * ON ENTRY :
 *   datare    real parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   dataim    imaginary parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the addition jobs;
 *   in2idx    indices of the second input of the addition jobs;
 *   outidx    indices of the output of the addition jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   datare    updated real parts of 
 *             the forward, backward, and cross products;
 *   dataim    updated imaginary parts of 
 *             the forward, backward, and cross products. */

void dbl_convoluted_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   data     coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products;
 *   output   space for the value and all derivatives;
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
 *   output   contains the value and all derivatives. */

void cmplx_convoluted_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datare   real parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   dataim   imaginary parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   outputre has space for all real parts of the output; 
 *   outputim has space for all imaginary parts of the output; 
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
 *   outputre stores the real parts of the value and all derivatives;
 *   outputim stores the imaginary parts of the value and all derivatives. */

void dbl_added_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   data     coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   output   space for the value and all derivatives;
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
 *   jobs     defines all addition jobs;
 *   verbose  if true, writes extra information.
 *
 * ON RETURN :
 *   output   contains the value and all derivatives. */

void cmplx_added_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr, int **idx,
   int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datare   real parts of coefficients of monomials and input series, 
 *            computed real parts of forward, backward, and cross products,
 *            with the accumulated additions;
 *   dataim   imaginary parts of coefficients of monomials and input series, 
 *            computed imaginary parts of forward, backward, and cross
 *            products, with the accumulated additions;
 *   outputre has space for all real parts of the output; 
 *   outputim has space for all imaginary parts of the output; 
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
 *   jobs     defines all addition jobs;
 *   verbose  if true, writes extra information.
 *
 * ON RETURN :
 *   outputre stores the real parts of the value and all derivatives;
 *   outputim stores the imaginary parts of the value and all derivatives. */

void write_coefficient_indices
 ( int cnt, int nbr, int *fsums, int *fstart, int *bsums, int *bstart,
   int *csums, int *cstart );
/*
 * DESCRIPTION :
 *   Writes the information about the coefficient indices in the data.
 *
 * ON ENTRY :
 *   cnt      total coefficient count,
 *   nbr      number of monomials, excluding the constant term;
 *   fsums    fsums[k] holds the sum of coefficients for the forward
 *            products of the k-th monomial and of all monomials before k;
 *   bsums    fsums[k] holds the sum of coefficients for the backward
 *            products of the k-th monomial and of all monomials before k;
 *   csums    fsums[k] holds the sum of coefficients for the cross
 *            products of the k-th monomial and of all monomials before k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial. */

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by jobs,
 *   on real data.
 *
 * ON ENTRY :
 *   BS       number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cst      deg+1 doubles for the constant coefficient series;
 *   cff      cff[k] has deg+1 doubles for the coefficient series
 *            of monomial k;
 *   input    contains the coefficients of the power series
 *            for all variables in the polynomial;
 *   output   space allocated for the value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   addjobs  addition jobs organized in layers;
 *   verbose  if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim, 
   double **inputre, double **inputim, double **outputre, double **outputim,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by jobs,
 *   on complex data.
 *
 * ON ENTRY :
 *   BS       number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cstre    deg+1 doubles for the real parts of the constant term;
 *   cstim    deg+1 doubles for the imaginary parts of the constant terms;
 *   cffre    cffre[k] has deg+1 doubles for the real parts of
 *            the coefficient series of monomial k;
 *   cffim    cffim[k] has deg+1 doubles for the imaginary parts of
 *            the coefficient series of monomial k;
 *   inputre  has the real parts of the coefficients of the power series
 *            for all variables in the polynomial;
 *   inputim  has the imaginary parts of the coefficients of the power series
 *            for all variables in the polynomial;
 *   outputre has space allocated for the real parts of
 *            the value and all derivatives;
 *   outputim has space allocated for the imaginary parts of
 *            the value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   addjobs  addition jobs organized in layers;
 *   verbose  if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputre stores the real parts of the derivatives and
 *            the value of the polynomial,
 *            outputre[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputre[dim] contains the value of the polynomial;
 *   outputim stores the imaginary parts of the derivatives and
 *            the of the polynomial,
 *            outputim[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputim[dim] contains the value of the polynomial;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

#endif
