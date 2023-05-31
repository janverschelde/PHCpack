// The file dbl_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in double precision.

#ifndef __dbl_polynomials_kernels_h__
#define __dbl_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"

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

__global__ void dbl_increment_jobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim );
/*
 * DESCRIPTION :
 *   Executes all increment jobs at the same layer.
 *   The block index defines the increment job.
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

__global__ void cmplxvectorized_flipsigns
 ( double *datari, int *flpidx, int dim );
/*
 * DESCRIPTION :
 *   Kernel to flip the signs of the second real operand
 *   on the data arrays used for the complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   datari       doubles of the convolutions;
 *   flpidx       start indices of series to flip;
 *   dim          equals the size of each block, or deg+1,
 *                where deg is the degree of truncation.
 *
 * ON RETURN :
 *   datari       doubles of the computed data. */

void GPU_cmplxvectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datari, double *elapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Flips the signs in the second operand of the real convolutions
 *   in the data arrays used in the vectorized complex arithmetic.
 *
 * ON ENTRY :
 *   deg          degree of truncation of the series;
 *   nbrflips     number of series to flip signs, number of blocks;
 *   flipidx      start index of every series to flip;
 *   datari       doubles of the convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datari       doubles of the computed data;
 *   elapsedms    elapsed time expressed in milliseconds. */

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
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
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
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   output   contains the value and all derivatives. */

void cmplx_convoluted_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
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
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   outputre stores the real parts of the value and all derivatives;
 *   outputim stores the imaginary parts of the value and all derivatives. */

void dbl_added_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   int vrblvl );
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
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   output   contains the value and all derivatives. */

void cmplx_added_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr, int **idx,
   int *fstart, int *bstart, int *cstart, AdditionJobs jobs, int vrblvl );
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
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   outputre stores the real parts of the value and all derivatives;
 *   outputim stores the imaginary parts of the value and all derivatives. */

void cmplx_added_datavectorized_to_output
 ( double *datari, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totcff, int offsetri, ComplexAdditionJobs jobs, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datarih      doubles of the computed data;
 *   outputre     has space for all doubles of the real parts
 *                of value and all derivatives;
 *   outputim     has space for all doubles of the imaginary parts
 *                of value and all derivatives;
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] is the number of variables for monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   totcff       total number of coefficients without vectorization;
 *   offsetri     size of the second operand;
 *   jobs         defines all addition jobs;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputre     contains the doubles of the real parts
 *                of the value and all derivatives;
 *   outputim     contains the doubles of the imaginary parts
 *                of the value and all derivatives. */

void dbl_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *data, double *cst, double **cff, double **input );
/*
 * DESCRIPTION :
 *   Initializes the real data vectors with the constants, the coefficients
 *   of the monomials and the input series for each variable.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   totcff       total number of coefficients;
 *   data         space for doubles of the data;
 *   cst          constant coefficient series;
 *   cff          cff[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   input        has the power series for all variables in the polynomial.
 *
 * ON RETURN :
 *   data         doubles of the initialized data. */

void cmplx_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datare, double *dataim, double *cstre, double *cstim,
   double **cffre, double **cffim, double **inputre, double **inputim );
/*
 * DESCRIPTION :
 *   Initializes the complex data vectors with the constants, the coefficients
 *   of the monomials and the input series for each variable.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   datare       space for doubles of the real data;
 *   dataim       space for doubles of the imag data;
 *   cstre        deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstim        deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffre        has the doubles of the real parts
 *                of the coefficients, cffrehi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffim        has the doubles of the imaginary parts
 *                of the coefficients, cffimhi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   inputre      has the doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputim      has the doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *
 * ON RETURN :
 *   datare       doubles of the initialized real data;
 *   dataim       doubles of the initialized imag data. */

void cmplxvectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datarihi, double *cstre, double *cstim,
   double **cffre, double **cffim, double **inputre, double **inputim );
/*
 * DESCRIPTION :
 *   Initializes the data vectors with the constants, the coefficients
 *   of the monomials and the input series for each variable,
 *   in data arrays where the imaginary parts follow the real parts.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   totcff       total number of coefficients without vectorization;
 *   offsetri     size of the second operand;
 *   datari       space for doubles of data;
 *   cstre        deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstim        deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffre        has the doubles of the real parts
 *                of the coefficients, cffrehi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffim        has the doubles of the imaginary parts
 *                of the coefficients, cffimhi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   inputre      has the doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputim      has the doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial.
 *
 * ON RETURN :
 *   datari       doubles of the initialized data. */

void dbl_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *data, double *cnvlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all convolution jobs on real data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   cnvjobs      convolution jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   data         doubles of the initialized data;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   data         doubles of the convolutions.
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds. */

void cmplx_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datare, double *dataim, double *cnvlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all convolution jobs on complex data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   cnvjobs      convolution jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   datare       doubles of the initialized real data;
 *   dataim       doubles of the initialized imag data;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datare       doubles of the real convolutions;
 *   dataim       doubles of the imag convolutions;
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds. */

void cmplxvectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datari, double *cnvlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all convolution jobs on complex data,
 *   on data arrays suitable for complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   totcff       total number of coefficients without vectorization;
 *   offsetri     size of the second operand;
 *   cnvjobs      convolution jobs organized in layers;
 *   incjobs      increment jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   datari       doubles of the initialized data;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datari       doubles of the convolutions;
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds. */

void dbl_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *data, double *addlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all addition jobs on real data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   addjobs      addition jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   data         doubles of the convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   data         doubles of the added convolutions;
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds. */

void cmplx_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datare, double *dataim, double *addlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all addition jobs on complex data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   addjobs      addition jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   datare       doubles of the real convolutions;
 *   dataim       doubles of the imag convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datare       doubles of the added real convolutions;
 *   dataim       doubles of the added imag convolutions;
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds. */

void cmplxvectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datari, double *addlapms, int vrblvl );
/*
 * DESCRIPTION :
 *   Launches the kernels for all addition jobs on complex data,
 *   on data arrays suitable for complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   totcff       total number of coefficients without vectorization;
 *   offsetri     size of the second operand;
 *   addjobs      addition jobs organized in layers;
 *   fstart       fstart[k] has the start position of the forward products
 *                for the k-th monomial;
 *   bstart       fstart[k] has the start position of the backward products
 *                for the k-th monomial;
 *   cstart       fstart[k] has the start position of the cross products
 *                for the k-th monomial;
 *   datari       doubles of the convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datari       doubles of the added convolutions.
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds. */

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl );
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
 *   vrblvl   is the verbose level.
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
   double *walltimesec, int vrblvl );
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
 *   vrblvl   is the verbose level.
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

void GPU_cmplxvectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim,
   double **inputre, double **inputim, double **outputre, double **outputim,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   ComplexAdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in several variables.
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   on complex data.  The complex arithmetic is vectorized,
 *   as defined by the complex versions of the jobs.
 *
 * ON ENTRY :
 *   BS             number of threads in a block, must equal deg + 1;
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   cstre          deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstim          deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cffre          has the doubles of the real parts
 *                  of the coefficients, cffrehihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   cffim          has the doubles of the imaginary parts
 *                  of the coefficients, cffimhihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   inputre        has the doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputim        has the doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   outputre       has space for the doubles of the real parts
 *                  of the value and all derivatives;
 *   outputim       has space for the doubles of the imaginary parts
 *                  of the value and all derivatives;
 *   cnvjobs        convolution jobs organized in layers;
 *   incjobs        increment jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputre       has the doubles of the real parts,
 *   outputim       has the doubles of the imaginary parts,
 *                  of derivatives and the value,
 *                  output[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  output[dim] contains the value of the polynomial;
 *   cnvlapms       is the elapsed time spent by all convolution kernels,
 *                  expressed in milliseconds;
 *   addlapms       is the elapsed time spent by all addition kernels,
 *                  expressed in milliseconds;
 *   elapsedms      is the elapsed time spent by all kernels,
 *                  expressed in milliseconds;
 *   walltimesec    is the elapsed wall clock time for all computations
 *                  (excluding memory copies) in seconds. */

#endif
