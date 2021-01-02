// The file dbl3_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in triple double precision.

#ifndef __dbl3_polynomials_kernels_h__
#define __dbl3_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "dbl3_polynomials_kernels.h"

__global__ void dbl3_padded_convjobs
 ( double *datahi, double *datami, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim );
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
 *   datahi   high parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datami   middle parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datalo   low parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   in1idx   indices of the first input of the convolution jobs;
 *   in2idx   indices of the second input of the convolution jobs;
 *   outidx   indices of the output of the convolution jobs;
 *   dim      the number of coefficients in each series
 *            equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahi   updated high forward, backward, and cross products;
 *   datami   updated middle forward, backward, and cross products;
 *   datalo   updated low forward, backward, and cross products. */

__global__ void dbl3_update_addjobs
 ( double *datahi, double *datami, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all addition jobs at the same layer.
 *   The block index defines the addition job.
 *
 * REQUIRED : 
 *   The number of blocks equals the size of  in1idx, in2idx, outidx,
 *   and dim equals the number of threads in each block.
 *
 * ON ENTRY :
 *   datahi   high parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datami   middle parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datalo   low parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   in1idx   indices of the first input of the addition jobs;
 *   in2idx   indices of the second input of the addition jobs;
 *   outidx   indices of the output of the addition jobs;
 *   dim      the number of coefficients in each series
 *            equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahi   updated high forward, backward, and cross products;
 *   datami   updated middle forward, backward, and cross products;
 *   datalo   updated low forward, backward, and cross products. */

void convoluted_data3_to_output
 ( double *datahi, double *datami, double *datalo,
   double **outputhi, double **outputmi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datahi   high parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datami   middle parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datalo   low parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   outputhi has space for the high parts of value and all derivatives;
 *   outputmi has space for the middle parts of value and all derivatives;
 *   outputlo has space for the low parts of value and all derivatives;
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
 *   outputhi has the high parts of the value and all derivatives;
 *   outputmi has the middle parts of the value and all derivatives;
 *   outputlo has the low parts of the value and all derivatives. */

void added_data3_to_output
 ( double *datahi, double *datami, double *datalo,
   double **outputhi, double **outputmi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datahi   high coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   datami   middle coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   datalo   low coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   outputhi has space for the high parts of value and all derivatives;
 *   outputmi has space for the middle parts of value and all derivatives;
 *   outputlo has space for the low parts of value and all derivatives;
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
 *   outputhi has the high parts of the value and all derivatives;
 *   outputmi has the middle parts of the value and all derivatives;
 *   outputlo has the low parts of the value and all derivatives. */

void GPU_dbl3_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi,  double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose=true );
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
 *   csthi    high deg+1 doubles for the constant coefficient series;
 *   cstmi    middle deg+1 doubles for the constant coefficient series;
 *   cstlo    low deg+1 doubles for the constant coefficient series;
 *   cffhi    cffhi[k] has deg+1 doubles for the high parts of the
 *            coefficient series of monomial k;
 *   cffmi    cffmi[k] has deg+1 doubles for the middle parts of the
 *            coefficient series of monomial k;
 *   cfflo    cfflo[k] has deg+1 doubles for the low parts of the
 *            coefficient series of monomial k;
 *   inputhi  has the high parts of coefficients of the power series
 *            for all variables in the polynomial;
 *   inputmi  has the middle parts of coefficients of the power series
 *            for all variables in the polynomial;
 *   inputlo  has the low parts of coefficients of the power series
 *            for all variables in the polynomial;
 *   outputhi has space for high parts of the value and all derivatives;
 *   outputmi has space for middle parts of the value and all derivatives;
 *   outputlo has space for low parts of the value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   addjobs  addition jobs organized in layers;
 *   verbose  if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputhi has the high parts of derivatives and the value,
 *            outputhi[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputhi[dim] contains the value of the polynomial;
 *   outputmi has the middle parts of derivatives and the value,
 *            outputmi[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputmi[dim] contains the value of the polynomial;
 *   outputlo has the low parts of derivatives and the value,
 *            outputhi[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputhi[dim] contains the value of the polynomial. */

#endif
