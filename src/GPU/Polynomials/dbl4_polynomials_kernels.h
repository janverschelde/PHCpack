// The file dbl4_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in quad double precision.

#ifndef __dbl4_polynomials_kernels_h__
#define __dbl4_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

__global__ void dbl4_padded_convjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
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
 *   datahihi  highest parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   datalohi  second highest parts of coefficients of monomials and input, 
 *             space for forward, backward, and cross products;
 *   datahilo  second lowest parts of coefficients of monomials and input, 
 *             space for forward, backward, and cross products;
 *   datalolo  lowest parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the convolution jobs;
 *   in2idx    indices of the second input of the convolution jobs;
 *   outidx    indices of the output of the convolution jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahihi  updated highest forward, backward, and cross products;
 *   datalohi  updated second highest forward, backward, and cross products;
 *   datahilo  updated second lowest forward, backward, and cross products;
 *   datalolo  updated lowest forward, backward, and cross products. */

__global__ void cmplx4_padded_convjobs
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
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
 *   datarehihi   highest doubles of the real parts of the coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelohi   second highest doubles of the real parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarehilo   second lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelolo   lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhihi   highest doubles of the imaginary parts 
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlohi   second highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhilo   second lowest doubles of the imaginary parts
 *                of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlolo   lowest doubles of the imaginary parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   in1idx       indices of the first input of the convolution jobs;
 *   in2idx       indices of the second input of the convolution jobs;
 *   outidx       indices of the output of the convolution jobs;
 *   dim          the number of coefficients in each series
 *                equals the number of threads in each block.
 *
 * ON RETURN :
 *   datarehihi   updated highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelohi   updated second highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarehilo   updated second lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelolo   updated lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   dataimhihi   updated highest doubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimlohi   updated second highest doubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimhilo   updated second lowest fdoubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimlolo   updated lowest fdoubles of the imaginary parts
 *                of the forward, backward, and cross products. */

__global__ void dbl4_update_addjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
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
 *   datahihi  highest parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   datalohi  second highest parts of coefficients of monomials and input, 
 *             space for forward, backward, and cross products;
 *   datahilo  second lowest parts of coefficients of monomials and input, 
 *             space for forward, backward, and cross products;
 *   datalolo  lowest parts of coefficients of monomials and input series, 
 *             space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the addition jobs;
 *   in2idx    indices of the second input of the addition jobs;
 *   outidx    indices of the output of the addition jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahihi  updated highest forward, backward, and cross products;
 *   datalohi  updated second highest forward, backward, and cross products;
 *   datahilo  updated second lowest forward, backward, and cross products;
 *   datalolo  updated lowest forward, backward, and cross products. */

__global__ void cmplx4_update_addjobs
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
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
 *   datarehihi   highest doubles of the real parts of the coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelohi   second highest doubles of the real parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarehilo   second lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelolo   lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhihi   highest doubles of the imaginary parts 
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlohi   second highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhilo   second lowest doubles of the imaginary parts
 *                of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlolo   lowest doubles of the imaginary parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   in1idx       indices of the first input of the addition jobs;
 *   in2idx       indices of the second input of the addition jobs;
 *   outidx       indices of the output of the addition jobs;
 *   dim          the number of coefficients in each series
 *                equals the number of threads in each block.
 *
 * ON RETURN :
 *   datarehihi   updated highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelohi   updated second highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarehilo   updated second lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelolo   updated lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   dataimhihi   updated highest doubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimlohi   updated second highest doubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimhilo   updated second lowest fdoubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimlolo   updated lowest fdoubles of the imaginary parts
 *                of the forward, backward, and cross products. */

void dbl_convoluted_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datahihi     highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datalohi     second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahilo     second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalolo     lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
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
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial. */

void cmplx_convoluted_data4_to_output
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datarehihi   highest doubles of the real parts of the coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelohi   second highest doubles of the real parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarehilo   second lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelolo   lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhihi   highest doubles of the imaginary parts 
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlohi   second highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhilo   second lowest doubles of the imaginary parts
 *                of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlolo   lowest doubles of the imaginary parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   outputrehihi has space for all highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelohi has space for all second highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrehilo has space for all second lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelolo has space for all lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputimhihi has space for all highest doubles of the imaginary parts
 *                of value and all derivatives;
 *   outputimlohi has space for all second highest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimhilo has space for all second lowest doubles of the
 *                imaginary parts of value and all derivatives;
 *   outputimlolo has space for all lowest doubles of the imaginary parts
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
 *   verbose      if true, writes extra information.
 *
 * ON RETURN :
 *   outpurethihi contains the highest doubles of the real parts
 *                of the value and all derivatives;
 *   outpuretlohi contains the second highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrehilo contains the second lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelolo contains the lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outpuimthihi contains the highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outpuimtlohi contains the second highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhilo contains the second lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlolo contains the lowest doubles of the imaginary parts
 *                of the value and all derivatives. */

void dbl_added_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datahihi     highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalohi     second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahilo     second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalolo     lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
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
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial. */

void cmplx_added_data4_to_output
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datarehihi   highest doubles of the real parts of the coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelohi   second highest doubles of the real parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarehilo   second lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datarelolo   lowest doubles of the real parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhihi   highest doubles of the imaginary parts 
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlohi   second highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimhilo   second lowest doubles of the imaginary parts
 *                of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   dataimlolo   lowest doubles of the imaginary parts of coefficients
 *                of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   outputrehihi has space for all highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelohi has space for all second highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrehilo has space for all second lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelolo has space for all lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputimhihi has space for all highest doubles of the imaginary parts
 *                of value and all derivatives;
 *   outputimlohi has space for all second highest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimhilo has space for all second lowest doubles of the
 *                imaginary parts of value and all derivatives;
 *   outputimlolo has space for all lowest doubles of the imaginary parts
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
 *   jobs         defines all addition jobs;
 *   verbose      if true, writes extra information.
 *
 * ON RETURN :
 *   outpurethihi contains the highest doubles of the real parts
 *                of the value and all derivatives;
 *   outpuretlohi contains the second highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrehilo contains the second lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelolo contains the lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outpuimthihi contains the highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outpuimtlohi contains the second highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhilo contains the second lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlolo contains the lowest doubles of the imaginary parts
 *                of the value and all derivatives. */

void GPU_dbl4_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose=true );
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
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second higest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihi    has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohi    has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilo    has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolo    has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   verbose      if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial;
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds;
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds;
 *   elapsedms    is the elapsed time spent by all kernels,
 *                expressed in milliseconds;
 *   walltimesec  is the elapsed wall clock time for all computations
 *                (excluding memory copies) in seconds. */

void GPU_cmplx4_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
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
 *   BS             number of threads in a block, must equal deg + 1;
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   cstrehihi      highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrelohi      second highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrehilo      second lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstrelolo      lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstimhihi      highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlohi      second highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimhilo      second lowest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlolo      lowest deg+1 doubles for the imaginary parts
 *                  of the constant coefficient series;
 *   cffrehihi      has the highest doubles of the real parts
 *                  of the coefficients, cffrehihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   cffrelohi      has the second highest doubles of the real parts
 *                  of the coefficients, cffrelohi[k] has deg+1 second highest
 *                  coefficients of monomial k;
 *   cffrehilo      has the second lowest doubles of the real parts
 *                  of the coefficients, cffrehilo[k] has deg+1 second lowest
 *                  coefficients of monomial k;
 *   cffrelolo      has the lowest doubles of the real parts
 *                  of the coefficients, cffrelolo[k] has deg+1 lowest
 *                  coefficients of monomial k;
 *   cffimhihi      has the highest doubles of the imaginary parts
 *                  of the coefficients, cffimhihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   cffimlohi      has the second highest doubles of the imaginary parts
 *                  of the coefficients, cffimlohi[k] has deg+1 second highest
 *                  coefficients of monomial k;
 *   cffimhilo      has the second lowest doubles of the imaginary parts
 *                  of the coefficient, cffimlolo[k] has the deg+1 second
 *                  lowest coefficients of monomial k;
 *   cffimlolo      has the lowest doubles of the imaginary parts
 *                  of the coefficient, cffimlolo[k] has the deg+1 lowest
 *                  coefficients of monomial k;
 *   inputrehihi    has the highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelohi    has the second highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelolo    has the second lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelolo    has the lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhihi    has the highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlohi    has the second highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhilo    has the second lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlolo    has the lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   outputrehihi   has space for the highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelohi   has space for the second highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrehilo   has space for the second lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelolo   has space for the lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputimhihi   has space for the highest doubles of the imaginary parts
 *                  of the value and all derivatives;
 *   outputimlohi   has space for the second highest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimhilo   has space for the second lowest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimlolo   has space for the lowest doubles of the imaginary parts
 *                  of the value and all derivatives;
 *   cnvjobs        convolution jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   verbose        if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputrehihi   has the highest doubles of the real parts,
 *   outputrelohi   has the second highest doubles of the real parts,
 *   outputrehilo   has the second lowest doubles of the real parts,
 *   outputrelolo   has the lowest doubles of the real parts,
 *   outputimhihi   has the highest doubles of the imaginary parts,
 *   outputimlohi   has the second highest doubles of the imaginary parts,
 *   outputimhilo   has the second lowest doubles of the imaginary parts,
 *   outputimlolo   has the lowest doubles of the imaginary parts
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
