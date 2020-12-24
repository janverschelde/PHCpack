// The file dbl_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in double precision.

#ifndef __dbl_polynomials_kernels_h__
#define __dbl_polynomials_kernels_h__

#include "convolution_jobs.h"

int coefficient_count ( int dim, int nbr, int deg, int *nvr );
/*
 * DESCRIPTION :
 *   Returns the total number of coefficients in all convolutions.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k. */

void coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart );
/*
 * DESCRIPTION :
 *   Computes the sums of coefficients for forward, backward, cross products,
 *   and the start positions for each monomial.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k.
 *
 * ON RETURN :
 *   fsums    fsums[k] holds the sum of coefficients for the forward
 *            products of the k-th monomial and of all monomials before k;
 *   bsums    fsums[k] holds the sum of coefficients for the backward
 *            products of the k-th monomial and of all monomials before k;
 *   csums    fsums[k] holds the sum of coefficients for the cross
 *            products of the k-th monomial and of all monomials before k;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial. */

void job_indices
 ( ConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Computes the indices of the two inputs and the output of a job.
 *
 * ON ENTRY :
 *   job      defines a convolution job;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the job.
 *
 * ON RETURN :
 *   inp1ix   index of the first input;
 *   inp2ix   index of the second input;
 *   outidx   index of the output. */

void jobs_coordinates
 ( ConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg,
   int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Defines the coordinates of all jobs in the same layer.
 *
 * ON ENTRY :
 *   jobs     defines convolution jobs;
 *   layer    the index of one layer of jobs;
 *   inp1ix   space for as many integers as the jobs on the layer;
 *   inp2ix   space for as many integers as the jobs on the layer;
 *   outidx   space for as many integers as the jobs on the layer;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   verbose  if true, writes extra information about the jobs.
 *
 * ON RETURN :
 *   inp1ix   inp1ix[i] is the index of the first input of job i;
 *   inp2ix   inp2ix[i] is the index of the second input of job i;
 *   outidx   outidx[i] is the index of the output of job i. */

__global__ void dbl_padded_convjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all convolution jobs at the same layer.
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

void data_to_output
 ( double *data, double *cst, double **output,
   int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   data     coefficients of monomials and input series, 
 *            computed forward, backward, and cross products;
 *   cst      constant coefficient of the polynomial;
 *   output   space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   fsums    fsums[k] holds the sum of coefficients for the forward
 *            products of the k-th monomial and of all monomials before k;
 *   bsums    fsums[k] holds the sum of coefficients for the backward
 *            products of the k-th monomial and of all monomials before k;
 *   csums    fsums[k] holds the sum of coefficients for the cross
 *            products of the k-th monomial and of all monomials before k;
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

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by jobs,
 *   all other parameters are the same as in the other function.
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
 *   jobs     convolution jobs organized in layers;
 *   verbose  if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

#endif
