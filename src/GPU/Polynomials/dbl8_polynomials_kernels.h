// The file dbl8_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in octo double precision.

#ifndef __dbl8_polynomials_kernels_h__
#define __dbl8_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"

__global__ void dbl8_padded_convjobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
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
 *   datalohihi  second highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilohi  third highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalolohi  fourth highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahihilo  fourth lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohilo  third lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilolo  second lowest parts of coefficients of monomials and input, 
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
 *   datalohihi  updated second highest forward, backward, and cross products;
 *   datahilohi  updated third highest forward, backward, and cross products;
 *   datalolohi  updated fourth highest forward, backward, and cross products;
 *   datahihilo  updated fourth lowest forward, backward, and cross products;
 *   datalohilo  updated third lowest forward, backward, and cross products;
 *   datahilolo  updated second lowest forward, backward, and cross products;
 *   datalololo  updated lowest forward, backward, and cross products. */

__global__ void dbl8_increment_jobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim );
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
 *   datahihihi are the highest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datalohihi are the 2nd highest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datahilohi are the 3rd highest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datalolohi are the 4th highest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   dataihhilo are the 4th lowest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datalohilo are the 3rd lowest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datahilolo are the 2nd lowest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   datalololo are the lowest parts of constant, coefficients,
 *             and input, space for forward, backward, and cross products;
 *   in1idx    indices of the first input of the increment jobs;
 *   in2idx    indices of the second input of the increment jobs;
 *   outidx    indices of the output of the increment jobs;
 *   dim       the number of coefficients in each series
 *             equals the number of threads in each block.
 *
 * ON RETURN :
 *   datahihihi has the incremented highest doubles;
 *   datalohihi has the incremented second highest doubles;
 *   datahilohi has the incremented third highest doubles;
 *   datalolohi has the incremented fourth highest doubles;
 *   datahihilo has the incremented fourth lowest doubles;
 *   datalohilo has the incremented third lowest doubles;
 *   datahilolo has the incremented second lowest doubles;
 *   datalololo has the incremented lowest doubles. */

__global__ void cmplx8_padded_convjobs
 ( double *datarehihihi, double *datarelohihi,
   double *datarehilohi, double *datarelolohi,
   double *datarehihilo, double *datarelohilo,
   double *datarehilolo, double *datarelololo,
   double *dataimhihihi, double *dataimlohihi,
   double *dataimhilohi, double *dataimlolohi,
   double *dataimhihilo, double *dataimlohilo,
   double *dataimhilolo, double *dataimlololo,
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
 *   datarehihihi are the highest doubles of the real parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarelohihi are the second highest doubles of the real parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarehilohi are the third highest doubles of the real parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarelolohi are the fourth highest doubles of the real parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarehihilo are the fourth lowest doubles of the real parts of
 *                the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarelohilo are the third lowest doubles of the real parts of
 *                the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarehilolo are the second lowest doubles of the real parts of
 *                the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   datarelololo are the lowest doubles of the real parts of 
 *                the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimhihihi are the highest doubles of the imaginary parts 
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimlohihi are the second highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimhilohi are the third highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimlolohi are the fourth highest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimhihilo are the fourth lowest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimlohilo are the third lowest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimhilolo are the second lowest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   dataimlololo are the lowest doubles of the imaginary parts
 *                of the coefficients of monomials and input series, with
 *                space for forward, backward, and cross products;
 *   in1idx       indices of the first input of the convolution jobs;
 *   in2idx       indices of the second input of the convolution jobs;
 *   outidx       indices of the output of the convolution jobs;
 *   dim          the number of coefficients in each series
 *                equals the number of threads in each block.
 *
 * ON RETURN :
 *   datarehihihi are the updated highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelohihi are the updated second highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarehilohi are the updated third highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelolohi are the updated fourth highest doubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarehihilo are the updated fourth lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelohilo are the updated third lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarehilolo are the updated second lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   datarelololo are the updated lowest fdoubles of the real parts
 *                of the forward, backward, and cross products;
 *   dataimhihihi are the updated highest doubles of the imaginary parts
 *                of the forward, backward, and cross products;
 *   dataimlohihi are the updated second highest doubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimhilohi are the updated third highest doubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimlolohi are the updated fourth highest doubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimhihilo are the updated fourth lowest fdoubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimlohilo are the updated third lowest fdoubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimhilolo are the updated second lowest fdoubles of the imaginary
 *                parts of the forward, backward, and cross products;
 *   dataimlololo are the updated lowest fdoubles of the imaginary parts
 *                of the forward, backward, and cross products. */

__global__ void cmplx8vectorized_flipsigns
 ( double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo, int *flpidx, int dim );
/*
 * DESCRIPTION :
 *   Kernel to flip the signs of the second real operand
 *   on the data arrays used for the complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   datarihihihi are the highest doubles of the convolutions;
 *   datarilohihi are the second highest doubles of the convolutions;
 *   datarihilohi are the third highest doubles of the convolutions;
 *   datarilolohi are the fourth highest doubles of the convolutions;
 *   datarihihilo are the fourth lowest doubles of the convolutions;
 *   datarilohilo are the third lowest doubles of the convolutions;
 *   datarihilolo are the second lowest doubles of the convolutions;
 *   datarilololo are the lowest doubles of the convolutions;
 *   flpidx       start indices of series to flip;
 *   dim          equals the size of each block, or deg+1,
 *                where deg is the degree of truncation.
 *
 * ON RETURN :
 *   datarihihihi are the highest doubles of the computed data;
 *   datarilohihi are the second highest doubles of the computed data;
 *   datarihilohi are the third highest doubles of the computed data;
 *   datarilolohi are the fourth highest doubles of the computed data;
 *   datarihihilo are the fourth lowest doubles of the computed data;
 *   datarilohilo are the third lowest doubles of the computed data;
 *   datarihilolo are the second lowest doubles of the computed data;
 *   datarilololo are the lowest doubles of the computed data. */

void GPU_cmplx8vectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double *elapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Flips the signs in the second operand of the real convolutions
 *   in the data arrays used in the vectorized complex arithmetic.
 *
 * ON ENTRY :
 *   deg          degree of truncation of the series;
 *   nbrflips     number of series to flip signs, number of blocks;
 *   flipidx      start index of every series to flip;
 *   datarihihihi are the highest doubles of the convolutions;
 *   datarilohihi are the second highest doubles of the convolutions;
 *   datarilohihi are the third highest doubles of the convolutions;
 *   datarilolohi are the fourth highest doubles of the convolutions;
 *   datarihihilo are the fourth lowest doubles of the convolutions;
 *   datarilohilo are the third lowest doubles of the convolutions;
 *   datarihilolo are the second lowest doubles of the convolutions;
 *   datarilololo are the lowest doubles of the convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datarihihihi are the highest doubles of the computed data;
 *   datarilohihi are the second highest doubles of the computed data;
 *   datarihilohi are the third highest doubles of the computed data;
 *   datarilolohi are the fourth highest doubles of the computed data;
 *   datarihihilo are the fourth lowest doubles of the computed data;
 *   datarilohilo are the third lowest doubles of the computed data;
 *   datarihilolo are the second lowest doubles of the computed data;
 *   datarilololo are the lowest doubles of the computed data;
 *   elapsedms    elapsed time expressed in milliseconds. */

__global__ void dbl8_update_addjobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
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
 *   datalohihi  second highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilohi  third highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalolohi  fourth highest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahihilo  fourth lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datalohilo  third lowest parts of coefficients of monomials and input, 
 *               space for forward, backward, and cross products;
 *   datahilolo  second lowest parts of coefficients of monomials and input, 
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
 *   datalohihi  updated second highest forward, backward, and cross products;
 *   datahilohi  updated third highest forward, backward, and cross products;
 *   datalolohi  updated fourth highest forward, backward, and cross products;
 *   datahihilo  updated fourth lowest forward, backward, and cross products;
 *   datalohilo  updated third lowest forward, backward, and cross products;
 *   datahilolo  updated second lowest forward, backward, and cross products;
 *   datalololo  updated lowest forward, backward, and cross products. */

void dbl_convoluted_data8_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datahihihi   highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   datalohihi   second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahilohi   third highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalolohi   fourth highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahihilo   fourth lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalohilo   third lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datahilolo   second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products;
 *   datalololo   lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlololo has space allocated for dim+1 series of degree deg;
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
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputlohihi has the second highest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputhilohi has the third highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputlolohi has the fourth highest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputhihilo has the fourth lowest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputlohilo has the third lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputhilolo has the second lowest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial. */

void dbl_added_data8_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datahihihi   highest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalohihi   second highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahilohi   third highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalolohi   fourth highest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahihilo   fourth lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalohilo   third lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datahilolo   second lowest parts of coefficients of monomials and input, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   datalololo   lowest parts of coefficients of monomials and input series, 
 *                space for forward, backward, and cross products,
 *                with the accumulated additions;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlololo has space allocated for dim+1 series of degree deg;
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
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputlohihi has the second highest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputhilohi has the third highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputlolohi has the fourth highest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputhihilo has the fourth lowest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputlohilo has the third lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputhilolo has the second lowest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial. */

void cmplx_added_data8vectorized_to_output
 ( double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totcff, int offsetri, ComplexAdditionJobs jobs, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datarihihihi has the highest doubles of the computed data;
 *   datarilohihi has the second highest doubles of the computed data;
 *   datarihilohi has the third highest doubles of the computed data;
 *   datarilolohi has the fourth highest doubles of the computed data;
 *   datarihihilo has the fourth lowest doubles of the computed data;
 *   datarilohilo has the third lowest doubles of the computed data;
 *   datarihilolo has the second lowest doubles of the computed data;
 *   datarilololo has the lowest doubles of the computed data;
 *   outputrehihihi has space for all highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelohihi has space for all second highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrehilohi has space for all third highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelolohi has space for all fourth highest doubles of the real parts
 *                of value and all derivatives;
 *   outputrehihilo has space for all fourth lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelohilo has space for all third lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputrehilolo has space for all second lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputrelololo has space for all lowest doubles of the real parts
 *                of value and all derivatives;
 *   outputimhihihi has space for all highest doubles of the imaginary parts
 *                of value and all derivatives;
 *   outputimlohihi has space for all second highest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimhilohi has space for all third highest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimlolohi has space for all fourth highest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimhihilo has space for all fourth lowest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimlohilo has space for all third lowest doubles of the 
 *                imaginary parts of value and all derivatives;
 *   outputimhilolo has space for all second lowest doubles of the
 *                imaginary parts of value and all derivatives;
 *   outputimlololo has space for all lowest doubles of the imaginary parts
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
 *   outputrehihihi has the highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelohihi has the second highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrehilohi has the third highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelolohi has the fourth highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrehihilo has the fourth lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelohilo has the third lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrehilolo has the second lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelololo has the lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputimhihihi has the highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlohihi has the second highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhilohi has the third highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhihihi has the fourth highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhihilo has the fourth lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlohilo has the third lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimhilolo has the second lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlololo has the lowest doubles of the imaginary parts
 *                of the value and all derivatives. */

void dbl8_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo );
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
 *   datahihihi   space for highest doubles of the data;
 *   datalohihi   space for the second highest doubles of the data;
 *   datahilohi   space for the third highest doubles of the data;
 *   datalolohi   space for the fourth highest doubles of the data;
 *   datahihilo   space for the fourth lowest doubles of the data;
 *   datalohilo   space for the third lowest doubles of the data;
 *   datahilolo   space for the second lowest doubles of the data;
 *   datalololo   space for the lowest doubles of the data;
 *   csthihihi    highest parts of constant coefficient series;
 *   cstlohihi    second highest parts of constant coefficient series;
 *   csthilohi    third highest parts of constant coefficient series;
 *   cstlolohi    fourth highest parts of constant coefficient series;
 *   csthihilo    fourth lowest parts of constant coefficient series;
 *   cstlohilo    third lowest parts of constant coefficient series;
 *   csthilolo    second lowest parts of constant coefficient series;
 *   cstlololo    lowest parts of constant coefficient series;
 *   cffhihihi    cffhihihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohihi    cfflohihi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilohi    cffhilohi[k] has deg+1 doubles for the third highest parts
 *                of the coefficient series of monomial k;
 *   cfflolohi    cfflolohi[k] has deg+1 doubles for the fourth highest parts
 *                of the coefficient series of monomial k;
 *   cffhihilo    cffhihilo[k] has deg+1 doubles for the fourth lowest parts
 *                of the coefficient series of monomial k;
 *   cfflohilo    cfflohilo[k] has deg+1 doubles for the third lowest parts
 *                of the coefficient series of monomial k;
 *   cffhilolo    cffhilolo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflololo    cfflololo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihihi  has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohihi  has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilohi  has the third highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolohi  has the fourth highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhihilo  has the fourth lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohilo  has the third lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilolo  has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlololo  has the lowest parts of the power series
 *                for all variables in the polynomial.
 *
 * ON RETURN :
 *   datahihihi   highest doubles of the initialized data;
 *   datalohihi   second highest doubles of the initialized data;
 *   datahilohi   third highest doubles of the initialized data;
 *   datalolohi   fourth highest doubles of the initialized data;
 *   datahihilo   fourth lowest doubles of the initialized data;
 *   datalohilo   third lowest doubles of the initialized data;
 *   datahilolo   second lowest doubles of the initialized data;
 *   datalololo   lowest doubles of the initialized data. */

void cmplx8vectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo );
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
 *   datarihihihi has space for highest doubles of data;
 *   datarilohihi has space for the second highest doubles of data;
 *   datarihilohi has space for the third highest doubles of data;
 *   datarilolohi has space for the fourth highest doubles of data;
 *   datarihihilo has space for the fourth lowest doubles of data;
 *   datarilohilo has space for the third lowest doubles of data;
 *   datarihilolo has space for the second lowest doubles of data;
 *   datarilololo has space for the lowest doubles of data;
 *   cstrehihihi  highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelohihi  second highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrehilohi  third highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelolohi  fourth highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrehihilo  fourth lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstrelohilo  third lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstrehilolo  second lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstrelololo  lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstimhihihi  highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlohihi  second highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimhilohi  third highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlolohi  fourth highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimhihilo  fourth lowest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlohilo  third lowest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimhilolo  second lowest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlololo  lowest deg+1 doubles for the imaginary parts
 *                of the constant coefficient series;
 *   cffrehihihi  has the highest doubles of the real parts
 *                of the coefficients, cffrehihihi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffrelohihi  has the second highest doubles of the real parts
 *                of the coefficients, cffrelohihi[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffrehilohi  has the third highest doubles of the real parts
 *                of the coefficients, cffrehilohi[k] has deg+1 third highest
 *                coefficients of monomial k;
 *   cffrelolohi  has the fourth highest doubles of the real parts
 *                of the coefficients, cffrelolohi[k] has deg+1 fourth highest
 *                coefficients of monomial k;
 *   cffrehihilo  has the fourth lowest doubles of the real parts
 *                of the coefficients, cffrehilolo[k] has deg+1 fourth lowest
 *                coefficients of monomial k;
 *   cffrelohilo  has the third lowest doubles of the real parts
 *                of the coefficients, cffrehilolo[k] has deg+1 third lowest
 *                coefficients of monomial k;
 *   cffrehilolo  has the second lowest doubles of the real parts
 *                of the coefficients, cffrehilolo[k] has deg+1 second lowest
 *                coefficients of monomial k;
 *   cffrelololo  has the lowest doubles of the real parts
 *                of the coefficients, cffrelololo[k] has deg+1 lowest
 *                coefficients of monomial k;
 *   cffimhihihi  has the highest doubles of the imaginary parts
 *                of the coefficients, cffimhihihi[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffimlohihi  has the second highest doubles of the imaginary parts
 *                of the coefficients, cffimlohihi[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimhilohi  has the third highest doubles of the imaginary parts
 *                of the coefficients, cffimhilohi[k] has deg+1 third highest
 *                coefficients of monomial k;
 *   cffimlolohi  has the fourth highest doubles of the imaginary parts
 *                of the coefficients, cffimlolohi[k] has deg+1 fourth highest
 *                coefficients of monomial k;
 *   cffimhilolo  has the third lowest doubles of the imaginary parts
 *                of the coefficient, cffimhilolo[k] has the deg+1 third
 *                lowest coefficients of monomial k;
 *   cffimlololo  has the lowest doubles of the imaginary parts
 *                of the coefficient, cffimlololo[k] has the deg+1 lowest
 *                coefficients of monomial k;
 *   inputrehihihi has the highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelohihi has the second highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrehilohi has the third highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelolohi has the fourth highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrehihilo has the fourth lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelohilo has the third lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrehilolo has the second lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelololo has the lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhihihi has the highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlohihi has the second highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhilohi has the third highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlolohi has the fourth highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhihilo has the fourth lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlohilo has the third lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhilolo has the second lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlololo has the lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial.
 *
 * ON RETURN :
 *   datarihihihi are the highest doubles of the initialized data;
 *   datarilohihi are the second highest doubles of the initialized data;
 *   datarihilohi are the third highest doubles of the initialized data;
 *   datarilolohi are the fourth highest doubles of the initialized data;
 *   datarihihilo are the fourth lowest doubles of the initialized data;
 *   datarilohilo are the third lowest doubles of the initialized data;
 *   datarihilolo are the second lowest doubles of the initialized data;
 *   datarilololo are the lowest doubles of the initialized data. */

void dbl8_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double *cnvlapms, int vrblvl );
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
 *   datahihihi   highest doubles of the initialized data;
 *   datalohihi   second highest doubles of the initialized data;
 *   datahilohi   third highest doubles of the initialized data;
 *   datalolohi   fourth highest doubles of the initialized data;
 *   datahihilo   fourth lowest doubles of the initialized data;
 *   datalohilo   third lowest doubles of the initialized data;
 *   datahilolo   second lowest doubles of the initialized data;
 *   datalololo   lowest doubles of the initialized data;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datahihihi   highest doubles of the convolutions;
 *   datalohihi   second highest doubles of the convolutions;
 *   datahilohi   third highest doubles of the convolutions;
 *   datalolohi   fourth highest doubles of the convolutions;
 *   datahihilo   fourth lowest doubles of the convolutions;
 *   datalohilo   third lowest doubles of the convolutions;
 *   datahilolo   second lowest doubles of the convolutions;
 *   datalololo   lowest doubles of the convolutions.
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds. */

void cmplx8vectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double *cnvlapms, int vrblvl );
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
 *   datarihihihi are the highest doubles of the initialized data;
 *   datarilohihi are the second highest doubles of the initialized data;
 *   datarihilohi are the third highest doubles of the initialized data;
 *   datarilolohi are the fourth highest doubles of the initialized data;
 *   datarihihilo are the fourth lowest doubles of the initialized data;
 *   datarilohilo are the third lowest doubles of the initialized data;
 *   datarihilolo are the second lowest doubles of the initialized data;
 *   datarilololo are the lowest doubles of the initialized data;
 *   vrblvl       ise the verbose level.
 *
 * ON RETURN :
 *   datarihihihi are the highest doubles of the convolutions;
 *   datarilohihi are the second highest doubles of the convolutions;
 *   datarihilohi are the third highest doubles of the convolutions;
 *   datarilolohi are the fourth highest doubles of the convolutions;
 *   datarihihilo are the fourth lowest doubles of the convolutions;
 *   datarilohilo are the third lowest doubles of the convolutions;
 *   datarihilolo are the second lowest doubles of the convolutions;
 *   datarilololo are the lowest doubles of the convolutions.
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds. */

void dbl8_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo, double *addlapms, int vrblvl );
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
 *   datahihihi   highest doubles of the convolutions;
 *   datalohihi   second highest doubles of the convolutions;
 *   datahilohi   third highest doubles of the convolutions;
 *   datalolohi   fourth highest doubles of the convolutions;
 *   datahihilo   fourth lowest doubles of the convolutions;
 *   datalohilo   third lowest doubles of the convolutions;
 *   datahilolo   second lowest doubles of the convolutions;
 *   datalololo   lowest doubles of the convolutions;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datahihihi   highest doubles of the added convolutions;
 *   datalohihi   second highest doubles of the added convolutions;
 *   datahilohi   third highest doubles of the added convolutions;
 *   datalolohi   fourth highest doubles of the added convolutions;
 *   datahihilo   fourth lowest doubles of the added convolutions;
 *   datalohilo   third lowest doubles of the added convolutions;
 *   datahilolo   second lowest doubles of the added convolutions;
 *   datalololo   lowest doubles of the added convolutions.
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds. */

void cmplx8vectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double *addlapms, int vrblvl );
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
 *   datarihihihi are the highest doubles of the convolutions;
 *   datarilohihi are the second highest doubles of the convolutions;
 *   datarihilohi are the third highest doubles of the convolutions;
 *   datarilolohi are the fourth highest doubles of the convolutions;
 *   datarihihilo are the fourth lowest doubles of the convolutions;
 *   datarilohilo are the third lowest doubles of the convolutions;
 *   datarihilolo are the second lowest doubles of the convolutions;
 *   datarilololo are the lowest doubles of the convolutions.
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   datarihihihi are the highest doubles of the added convolutions;
 *   datarilohihi are the second highest doubles of the added convolutions;
 *   datarihilohi are the third highest doubles of the added convolutions;
 *   datarilolohi are the fourth highest doubles of the added convolutions;
 *   datarihihilo are the fourth lowest doubles of the added convolutions;
 *   datarilohilo are the third lowest doubles of the added convolutions;
 *   datarihilolo are the second lowest doubles of the added convolutions;
 *   datarilololo are the lowest doubles of the added convolutions.
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds. */

void GPU_dbl8_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl );
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
 *   cstlohihi    second highest parts of constant coefficient series;
 *   csthilohi    third highest parts of constant coefficient series;
 *   cstlolohi    fourth highest parts of constant coefficient series;
 *   csthihilo    fourth lowest parts of constant coefficient series;
 *   cstlohilo    third lowest parts of constant coefficient series;
 *   csthilolo    second lowest parts of constant coefficient series;
 *   cstlololo    lowest parts of constant coefficient series;
 *   cffhihihi    cffhihihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohihi    cfflohihi[k] has deg+1 doubles for the second highest
 *                parts of the coefficient series of monomial k;
 *   cffhilohi    cffhilohi[k] has deg+1 doubles for the third highest
 *                parts of the coefficient series of monomial k;
 *   cfflolohi    cfflolohi[k] has deg+1 doubles for the fourth highest
 *                parts of the coefficient series of monomial k;
 *   cffhihilo    cffhihilo[k] has deg+1 doubles for the fourth lowest
 *                parts of the coefficient series of monomial k;
 *   cfflohilo    cfflohilo[k] has deg+1 doubles for the third lowest
 *                parts of the coefficient series of monomial k;
 *   cffhilolo    cffhilolo[k] has deg+1 doubles for the second lowest
 *                parts of the coefficient series of monomial k;
 *   cfflololo    cfflololo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihihi  has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohihi  has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilohi  has the third highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolohi  has the fourth highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhihilo  has the fourth lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohilo  has the third lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilolo  has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlololo  has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihihi has space allocated for dim+1 series of degree deg;
 *   outputlohihi has space allocated for dim+1 series of degree deg;
 *   outputhilohi has space allocated for dim+1 series of degree deg;
 *   outputlolohi has space allocated for dim+1 series of degree deg;
 *   outputhihilo has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
 *   outputhilolo has space allocated for dim+1 series of degree deg;
 *   outputlohilo has space allocated for dim+1 series of degree deg;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputhihihi has the highest parts of derivatives and the value,
 *                outputhihihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihihi[dim] contains the value of the polynomial;
 *   outputlohihi has the second highest parts of derivatives and the value,
 *                outputlohihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohihi[dim] contains the value of the polynomial;
 *   outputhilohi has the third highest parts of derivatives and the value,
 *                outputhilohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilohi[dim] contains the value of the polynomial;
 *   outputlolohi has the fourth highest parts of derivatives and the value,
 *                outputlolohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolohi[dim] contains the value of the polynomial;
 *   outputhihilo has the fourth lowest parts of derivatives and the value,
 *                outputhihilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihilo[dim] contains the value of the polynomial;
 *   outputlohilo has the third lowest parts of derivatives and the value,
 *                outputlohilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohilo[dim] contains the value of the polynomial;
 *   outputhilolo has the second lowest parts of derivatives and the value,
 *                outputhilolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilolo[dim] contains the value of the polynomial;
 *   outputlololo has the lowest parts of derivatives and the value,
 *                outputlololo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlololo[dim] contains the value of the polynomial;
 *   cnvlapms     is the elapsed time spent by all convolution kernels,
 *                expressed in milliseconds;
 *   addlapms     is the elapsed time spent by all addition kernels,
 *                expressed in milliseconds;
 *   elapsedms    is the elapsed time spent by all kernels,
 *                expressed in milliseconds;
 *   walltimesec  is the elapsed wall clock time for all computations
 *                (excluding memory copies) in seconds. */

void GPU_cmplx8vectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo,
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
 *   cstrehihihi    highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrelohihi    second highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrehilohi    third highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrelolohi    fourth highest deg+1 doubles of the real parts
 *                  of the constant coefficient series;
 *   cstrehihilo    fourth lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstrelohilo    third lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstrehilolo    second lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstrelololo    lowest deg+1 doubles for the real parts
 *                  of the constant coefficient series;
 *   cstimhihihi    highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlohihi    second highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimhilohi    third highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlolohi    fourth highest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimhihilo    fourth lowest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlohilo    third lowest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimhilolo    second lowest deg+1 doubles of the imaginary parts
 *                  of the constant coefficient series;
 *   cstimlololo    lowest deg+1 doubles for the imaginary parts
 *                  of the constant coefficient series;
 *   cffrehihihi    has the highest doubles of the real parts
 *                  of the coefficients, cffrehihihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   cffrelohihi    has the second highest doubles of the real parts
 *                  of the coefficients, cffrelohihi[k] has deg+1 second
 *                  highest coefficients of monomial k;
 *   cffrehilohi    has the third highest doubles of the real parts
 *                  of the coefficients, cffrehilohi[k] has deg+1 third
 *                  highest coefficients of monomial k;
 *   cffrelolohi    has the fourth highest doubles of the real parts
 *                  of the coefficients, cffrelolohi[k] has deg+1 fourth
 *                  highest coefficients of monomial k;
 *   cffrehihilo    has the fourth lowest doubles of the real parts
 *                  of the coefficients, cffrehihilo[k] has deg+1 fourth
 *                  lowest coefficients of monomial k;
 *   cffrelohilo    has the third lowest doubles of the real parts
 *                  of the coefficients, cffrelohilo[k] has deg+1 third
 *                  lowest coefficients of monomial k;
 *   cffrehilolo    has the second lowest doubles of the real parts
 *                  of the coefficients, cffrehilolo[k] has deg+1 second
 *                  lowest coefficients of monomial k;
 *   cffrelololo    has the lowest doubles of the real parts
 *                  of the coefficients, cffrelololo[k] has deg+1 lowest
 *                  coefficients of monomial k;
 *   cffimhihihi    has the highest doubles of the imaginary parts
 *                  of the coefficients, cffimhihihi[k] has deg+1 highest
 *                  coefficients of monomial k;
 *   cffimlohihi    has the second highest doubles of the imaginary parts
 *                  of the coefficients, cffimlohihi[k] has deg+1 second
 *                  highest coefficients of monomial k;
 *   cffimhilohi    has the third highest doubles of the imaginary parts
 *                  of the coefficients, cffimhilohi[k] has deg+1 third
 *                  highest coefficients of monomial k;
 *   cffimlolohi    has the fourth highest doubles of the imaginary parts
 *                  of the coefficients, cffimlolohi[k] has deg+1 fourth
 *                  highest coefficients of monomial k;
 *   cffimhihilo    has the fourth lowest doubles of the imaginary parts
 *                  of the coefficient, cffimhihilo[k] has the deg+1 fourth
 *                  lowest coefficients of monomial k;
 *   cffimlohilo    has the third lowest doubles of the imaginary parts
 *                  of the coefficient, cffimlohilo[k] has the deg+1 third
 *                  lowest coefficients of monomial k;
 *   cffimhilolo    has the second lowest doubles of the imaginary parts
 *                  of the coefficient, cffimhilolo[k] has the deg+1 second
 *                  lowest coefficients of monomial k;
 *   cffimlololo    has the lowest doubles of the imaginary parts
 *                  of the coefficient, cffimlololo[k] has the deg+1 lowest
 *                  coefficients of monomial k;
 *   inputrehihihi  has the highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelohihi  has the second highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrehilohi  has the third highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelolohi  has the fourth highest doubles of the real parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrehihilo  has the fourth lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelohilo  has the third lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrehilolo  has the second lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputrelololo  has the lowest doubles of the real part
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhihihi  has the highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlohihi  has the second highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhilohi  has the third highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlolohi  has the fourth highest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhihilo  has the fourth lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlohilo  has the third lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimhilolo  has the second lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   inputimlololo  has the lowest doubles of the imaginary parts
 *                  of the coefficients of the power series
 *                  for all variables in the polynomial;
 *   outputrehihihi has space for the highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelohihi has space for the second highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrehilohi has space for the third highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelolohi has space for the fourth highest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrehihilo has space for the fourth lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelohilo has space for the third lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrehilolo has space for the second lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputrelololo has space for the lowest doubles of the real parts
 *                  of the value and all derivatives;
 *   outputimhihihi has space for the highest doubles of the imaginary parts
 *                  of the value and all derivatives;
 *   outputimlohihi has space for the second highest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimhilohi has space for the third highest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimlolohi has space for the fourth highest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimhihilo has space for the fourth lowest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimlohilo has space for the third lowest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimhilolo has space for the second lowest doubles of
 *                  the imaginary parts of the value and all derivatives;
 *   outputimlololo has space for the lowest doubles of the imaginary parts
 *                  of the value and all derivatives;
 *   cnvjobs        convolution jobs organized in layers;
 *   incjobs        increment jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputrehihihi has the highest doubles of the real parts,
 *   outputrelohihi has the second highest doubles of the real parts,
 *   outputrehilohi has the third highest doubles of the real parts,
 *   outputrelolohi has the fourth highest doubles of the real parts,
 *   outputrehihilo has the fourth lowest doubles of the real parts,
 *   outputrelohilo has the third lowest doubles of the real parts,
 *   outputrehilolo has the second lowest doubles of the real parts,
 *   outputrelololo has the lowest doubles of the real parts,
 *   outputimhihihi has the highest doubles of the imaginary parts,
 *   outputimlohihi has the second highest doubles of the imaginary parts,
 *   outputimhilohi has the third highest doubles of the imaginary parts,
 *   outputimlolohi has the fourth highest doubles of the imaginary parts,
 *   outputimhihilo has the fourth lowest doubles of the imaginary parts,
 *   outputimlohilo has the third lowest doubles of the imaginary parts,
 *   outputimhilolo has the second lowest doubles of the imaginary parts,
 *   outputimlololo has the lowest doubles of the imaginary parts
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
