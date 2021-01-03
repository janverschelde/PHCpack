// The file dbl8_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in octo double precision.

#ifndef __dbl8_polynomials_kernels_h__
#define __dbl8_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "dbl8_polynomials_kernels.h"

__global__ void dbl8_padded_convjobs
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
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
 *   datahihihi  highest parts of coefficients of monomials and input series, 
 *               space for forward, backward, and cross products;
 *   datahilohi  second highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahihilo  third highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilolo  fourth highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohihi  fourth lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalolohi  third lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohilo  second lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalololo  lowest parts of coefficients of monomials and input series, 
 *               space for forward, backward, and cross products;
 *   in1idx      indices of the first input of the convolution jobs;
 *   in2idx      indices of the second input of the convolution jobs;
 *   outidx      indices of the output of the convolution jobs;
 *   dim         the number of coefficients in each series
 *               equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahihihi  updated highest forward, backward, and cross products;
 *   datahilohi  updated second highest forward, backward, and cross products;
 *   datahihilo  updated third highest forward, backward, and cross products;
 *   datahilolo  updated fourth highest forward, backward, and cross products;
 *   datalohihi  updated fourth lowest forward, backward, and cross products;
 *   datalolohi  updated third lowest forward, backward, and cross products;
 *   datalohilo  updated second lowest forward, backward, and cross products;
 *   datalololo  updated lowest forward, backward, and cross products. */

__global__ void dbl8_update_addjobs
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
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
 *   datahihihi  highest parts of coefficients of monomials and input series, 
 *               space for forward, backward, and cross products;
 *   datahilohi  second highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahihilo  third highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilolo  fourth highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohihi  fourth lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalolohi  third lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohilo  second lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalololo  lowest parts of coefficients of monomials and input series, 
 *               space for forward, backward, and cross products;
 *   in1idx      indices of the first input of the addition jobs;
 *   in2idx      indices of the second input of the addition jobs;
 *   outidx      indices of the output of the addition jobs;
 *   dim         the number of coefficients in each series
 *               equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahihihi  updated highest forward, backward, and cross products;
 *   datahilohi  updated second highest forward, backward, and cross products;
 *   datahihilo  updated third highest forward, backward, and cross products;
 *   datahilolo  updated fourth highest forward, backward, and cross products;
 *   datalohihi  updated fourth lowest forward, backward, and cross products;
 *   datalolohi  updated third lowest forward, backward, and cross products;
 *   datalohilo  updated second lowest forward, backward, and cross products;
 *   datalololo  updated lowest forward, backward, and cross products. */

void convoluted_data8_to_output
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datahihihi   highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datahilohi   second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahihilo   third highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahilolo   fourth highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalohihi   fourth lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalolohi   third lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalohilo   second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalololo   lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
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
 *   verbose      if true, writes extra information.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi has the second highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo has the third highest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo has the fourth highest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi has the fourth lowest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi has the third lowest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo has the second lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial. */

void added_data8_to_output
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datahihihi   highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahilohi   second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahihilo   third highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahilolo   fourth highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalohihi   fourth lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalolohi   third lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalohilo   second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalololo   lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
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
 *   jobs         defines all addition jobs;
 *   verbose      if true, writes extra information.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi has the second highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo has the third highest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo has the fourth highest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi has the fourth lowest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi has the third lowest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo has the second lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial. */

void GPU_dbl8_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo,
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Evaluates and differentiations a polynomial in 
 *   Computes the convolutions in the order as defined by jobs,
 *   all other parameters are the same as in the other function.
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
 *   csthihihi    highest parts of constant coefficient series;
 *   csthilohi    second higest parts of constant coefficient series;
 *   csthihilo    third higest parts of constant coefficient series;
 *   csthilolo    fourth higest parts of constant coefficient series;
 *   cstlohihi    fourth lowest parts of constant coefficient series;
 *   cstlolohi    third lowest parts of constant coefficient series;
 *   cstlohilo    second lowest parts of constant coefficient series;
 *   cstlololo    lowest parts of constant coefficient series;
 *   cffhihihi    cffhihihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cffhilohi    cffhilohi[k] has deg+1 doubles for the second highest
 *                parts of the coefficient series of monomial k;
 *   cffhihilo    cffhihilo[k] has deg+1 doubles for the third highest
 *                parts of the coefficient series of monomial k;
 *   cffhilolo    cffhilolo[k] has deg+1 doubles for the fourth highest
 *                parts of the coefficient series of monomial k;
 *   cfflohihi    cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                parts of the coefficient series of monomial k;
 *   cfflolohi    cfflolohi[k] has deg+1 doubles for the third lowest
 *                parts of the coefficient series of monomial k;
 *   cfflohilo    cfflohilo[k] has deg+1 doubles for the second lowest
 *                parts of the coefficient series of monomial k;
 *   cfflololo    cfflololo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihihi  has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilohi  has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhihilo  has the third highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilolo  has the fourth highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohihi  has the fourth lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolohi  has the third lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohilo  has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlololo  has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   verbose      if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi has the second highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo has the third highest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo has the fourth highest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi has the fourth lowest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi has the third lowest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo has the second lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial. */

#endif
