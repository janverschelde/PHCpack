// The file dbl5_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in penta double precision.

#ifndef __dbl5_polynomials_kernels_h__
#define __dbl5_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

__global__ void dbl5_padded_convjobs
 ( double *datatb, double *dataix, double *datami, 
   double *datarg, double *datapk,
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
 *   datatb  highest parts of coefficients of monomials and input series, 
 *           space for forward, backward, and cross products;
 *   dataix  second highest parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datami  middle parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datarg  second lowest parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datapk  lowest parts of coefficients of monomials and input series, 
 *           space for forward, backward, and cross products;
 *   in1idx  indices of the first input of the convolution jobs;
 *   in2idx  indices of the second input of the convolution jobs;
 *   outidx  indices of the output of the convolution jobs;
 *   dim     the number of coefficients in each series
 *           equals the number of threads in each block.
 *
 * ON RETURN :
 *   datatb  updated highest forward, backward, and cross products;
 *   dataix  updated second highest forward, backward, and cross products;
 *   datami  updated middle forward, backward, and cross products;
 *   datarg  updated second lowest forward, backward, and cross products;
 *   datapk  updated lowest forward, backward, and cross products. */

__global__ void dbl5_update_addjobs
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
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
 *   datatb  highest parts of coefficients of monomials and input series, 
 *           space for forward, backward, and cross products;
 *   dataix  second highest parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datami  middle parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datarg  second lowest parts of coefficients of monomials and input, 
 *           space for forward, backward, and cross products;
 *   datapk  lowest parts of coefficients of monomials and input series, 
 *           space for forward, backward, and cross products;
 *   in1idx  indices of the first input of the addition jobs;
 *   in2idx  indices of the second input of the addition jobs;
 *   outidx  indices of the output of the addition jobs;
 *   dim     the number of coefficients in each series
 *           equals the number of threads in each block.
 *
 * ON RETURN :
 *   datatb  updated highest forward, backward, and cross products;
 *   dataix  updated second highest forward, backward, and cross products;
 *   datami  updated middle forward, backward, and cross products;
 *   datarg  updated second lowest forward, backward, and cross products;
 *   datapk  updated lowest forward, backward, and cross products. */

void convoluted_data5_to_output
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datatb     highest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   dataix     second highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datami     middle parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datarg     second lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datapk     lowest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
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
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial. */

void added_data5_to_output
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datatb     highest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   dataix     second highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datami     middle parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datarg     second lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datapk     lowest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
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
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial. */

void GPU_dbl5_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
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
 *   BS         number of threads in a block, must equal deg + 1;
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second higest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   inputtb    has the highest parts of the power series
 *              for all variables in the polynomial;
 *   inputix    has the second highest parts of the power series
 *              for all variables in the polynomial;
 *   inputmi    has the middle parts of the power series
 *              for all variables in the polynomial;
 *   inputrg    has the second lowest parts of the power series
 *              for all variables in the polynomial;
 *   inputpk    has the lowest parts of the power series
 *              for all variables in the polynomial;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
 *   cnvjobs    convolution jobs organized in layers;
 *   addjobs    addition jobs organized in layers;
 *   verbose    if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial;
 *   cnvlapms   is the elapsed time spent by all convolution kernels,
 *              expressed in milliseconds;
 *   addlapms   is the elapsed time spent by all addition kernels,
 *              expressed in milliseconds;
 *   elapsedms  is the elapsed time spent by all kernels,
 *              expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *              (excluding memory copies) in seconds. */

#endif
