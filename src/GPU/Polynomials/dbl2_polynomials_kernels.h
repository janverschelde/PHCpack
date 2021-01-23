// The file dbl2_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in double double precision.

#ifndef __dbl2_polynomials_kernels_h__
#define __dbl2_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

__global__ void dbl2_padded_convjobs
 ( double *datahi, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim );
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
 *   datahi   high parts of coefficients of monomials and input series, 
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
 *   datalo   updated low forward, backward, and cross products. */

__global__ void cmplx2_padded_convjobs
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
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
 *   datarehi   high doubles of the real parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   datarelo   low doubles of the real parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimhi   high doubles of the imaginary parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimlo   low doubles of the imaginary parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   in1idx     indices of the first input of the convolution jobs;
 *   in2idx     indices of the second input of the convolution jobs;
 *   outidx     indices of the output of the convolution jobs;
 *   dim        the number of coefficients in each series
 *              equals the number of threads in each block.
 *
 * ON RETURN :
 *   datarehi   updated high doubles of the real parts
 *              of the forward, backward, and cross products;
 *   datarelo   updated low fdoubles of the real parts
 *              of the forward, backward, and cross products;
 *   dataimhi   updated high doubles of the imaginary parts
 *              of the forward, backward, and cross products;
 *   dataimlo   updated low fdoubles of the imaginary parts
 *              of the forward, backward, and cross products. */

__global__ void dbl2_update_addjobs
 ( double *datahi, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim );
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
 *   datahi   high parts of coefficients of monomials and input series, 
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
 *   datalo   updated low forward, backward, and cross products. */

__global__ void cmplx2_update_addjobs
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
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
 *   datarehi   high doubles of the real parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   datarelo   low doubles of the real parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimhi   high doubles of the imaginary parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimlo   low doubles of the imaginary parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   in1idx     indices of the first input of the addition jobs;
 *   in2idx     indices of the second input of the addition jobs;
 *   outidx     indices of the output of the addition jobs;
 *   dim        the number of coefficients in each series
 *              equals the number of threads in each block.
 *
 * ON RETURN :
 *   datarehi   updated high doubles of the real parts
 *              of the forward, backward, and cross products;
 *   datarelo   updated low fdoubles of the real parts
 *              of the forward, backward, and cross products;
 *   dataimhi   updated high doubles of the imaginary parts
 *              of the forward, backward, and cross products;
 *   dataimlo   updated low fdoubles of the imaginary parts
 *              of the forward, backward, and cross products. */

void dbl_convoluted_data2_to_output
 ( double *datahi, double *datalo, double **outputhi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datahi   high parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datalo   low parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   outputhi has space for the high parts of value and all derivatives;
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
 *   outputhi contains the high doubles of the value and all derivatives;
 *   outputlo contains the low doubles of the value and all derivatives. */

void cmplx_convoluted_data2_to_output
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datarehi   high doubles of the real parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   datarelo   low doubles of the real parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimhi   high doubles of the imaginary parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimlo   low doubles of the imaginary parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   outputrehi has space for all high doubles of the real parts
 *              of value and all derivatives;
 *   outputrelo has space for all low doubles of the real parts
 *              of value and all derivatives;
 *   outputimhi has space for all high doubles of the imaginary parts
 *              of value and all derivatives;
 *   outputimlo has space for all low doubles of the imaginary parts
 *              of value and all derivatives;
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] is the number of variables for monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   fstart     fstart[k] has the start position of the forward products
 *              for the k-th monomial;
 *   bstart     fstart[k] has the start position of the backward products
 *              for the k-th monomial;
 *   cstart     fstart[k] has the start position of the cross products
 *              for the k-th monomial;
 *   verbose    if true, writes extra information.
 *
 * ON RETURN :
 *   outpurethi contains the high doubles of the real parts
 *              of the value and all derivatives;
 *   outputrelo contains the low doubles of the real parts
 *              of the value and all derivatives;
 *   outpuimthi contains the high doubles of the imaginary parts
 *              of the value and all derivatives;
 *   outputimlo contains the low doubles of the imaginary parts
 *              of the value and all derivatives. */

void dbl_added_data2_to_output
 ( double *datahi, double *datalo, double **outputhi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datahi   high coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   datalo   low coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products,
 *            with the accumulated additions;
 *   outputhi has space for the high parts of value and all derivatives;
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
 *   outputlo has the low parts of the value and all derivatives. */

void cmplx_added_data2_to_output
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datarehi   high doubles of the real parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   datarelo   low doubles of the real parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimhi   high doubles of the imaginary parts of the coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataimlo   low doubles of the imaginary parts of coefficients
 *              of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   outputrehi has space for all high doubles of the real parts
 *              of value and all derivatives;
 *   outputrelo has space for all low doubles of the real parts
 *              of value and all derivatives;
 *   outputimhi has space for all high doubles of the imaginary parts
 *              of value and all derivatives;
 *   outputimlo has space for all low doubles of the imaginary parts
 *              of value and all derivatives;
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] is the number of variables for monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   fstart     fstart[k] has the start position of the forward products
 *              for the k-th monomial;
 *   bstart     fstart[k] has the start position of the backward products
 *              for the k-th monomial;
 *   cstart     fstart[k] has the start position of the cross products
 *              for the k-th monomial;
 *   jobs       defines all addition jobs;
 *   verbose    if true, writes extra information.
 *
 * ON RETURN :
 *   outpurethi contains the high doubles of the real parts
 *              of the value and all derivatives;
 *   outputrelo contains the low doubles of the real parts
 *              of the value and all derivatives;
 *   outpuimthi contains the high doubles of the imaginary parts
 *              of the value and all derivatives;
 *   outputimlo contains the low doubles of the imaginary parts
 *              of the value and all derivatives. */

void GPU_dbl2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
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
 *   csthi    high deg+1 doubles for the constant coefficient series;
 *   cstlo    low deg+1 doubles for the constant coefficient series;
 *   cffhi    cffhi[k] has deg+1 doubles for the high parts of the
 *            coefficient series of monomial k;
 *   cfflo    cfflo[k] has deg+1 doubles for the low parts of the
 *            coefficient series of monomial k;
 *   inputhi  has the high parts of coefficients of the power series
 *            for all variables in the polynomial;
 *   inputlo  has the low parts of coefficients of the power series
 *            for all variables in the polynomial;
 *   outputhi has space for high parts of the value and all derivatives;
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
 *   outputlo has the low parts of derivatives and the value,
 *            outputhi[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            outputhi[dim] contains the value of the polynomial.
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by jobs,
 *   performs the updates to the values as defined by addjobs,
 *   on complex data.
 *
 * ON ENTRY :
 *   BS           number of threads in a block, must equal deg + 1;
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cstrehi      high doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low doubles of the real parts
 *                of the constant coefficient series;
 *   cstimhi      high doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffrehi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffrelo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   inputrehi    has the high doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelo    has the low doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhi    has the high doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlo    has the low doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputrehi   has space for the high doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelo   has space for the low doubles of the real parts
 *                of the value and all derivatives;
 *   outputimhi   has space for the high doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlo   has space for the low doubles of the imaginary parts
 *                of the value and all derivatives;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   verbose      if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputrehi   has the high doubles of the real parts
 *                of the value and the derivatives,
 *   outputrelo   has the low doubles of the real parts
 *                of the value and the derivatives,
 *   outputimhi   has the high doubles of the imaginary parts
 *                of the value and the derivatives,
 *   outputimlo   has the low doubles of the imaginary parts
 *                of the value and the derivatives,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial;
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds;
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds;
 *   elapsedms    is the elapsed time spent by all kernels,
 *                expressed in milliseconds;
 *   walltimesec  is the elapsed wall clock time for all computations
 *                (excluding memory copies) in seconds. */

#endif
