// The file dbl10_ponomials_kernels.h specifies functions to evaluate and
// differentiate a ponomial at power series truncated to the same degree,
// in deca double precision.

#ifndef __dbl10_polynomials_kernels_h__
#define __dbl10_polynomials_kernels_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

__global__ void dbl10_padded_convjobs
 ( double *datartb, double *datarix, double *datarmi, 
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi, 
   double *datalrg, double *datalpk,
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
 *   datartb  highest parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datarix  second highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarmi  third highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarrg  fourth highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarpk  fifth highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   dataltb  fifth lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalix  fourth lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalmi  third lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalrg  second lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalpk  lowest parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   in1idx   indices of the first input of the convolution jobs;
 *   in2idx   indices of the second input of the convolution jobs;
 *   outidx   indices of the output of the convolution jobs;
 *   dim      the number of coefficients in each series
 *            equals the number of threads in each block.
 *
 * ON RETURN :
 *   datartb  updated highest forward, backward, and cross products;
 *   datarix  updated second highest forward, backward, and cross products;
 *   datarmi  updated third highest forward, backward, and cross products;
 *   datarrg  updated fourth highest forward, backward, and cross products;
 *   datarpk  updated fifth highest forward, backward, and cross products;
 *   dataltb  updated fifth lowest forward, backward, and cross products;
 *   datalix  updated fourth lowest forward, backward, and cross products;
 *   datalmi  updated third lowest forward, backward, and cross products;
 *   datalrg  updated second lowest forward, backward, and cross products;
 *   datalpk  updated lowest forward, backward, and cross products. */

__global__ void dbl10_update_addjobs
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
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
 *   datartb  highest parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   datarix  second highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarmi  third highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarrg  fourth highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datarpk  fifth highest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   dataltb  fifth lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalix  fourth lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalmi  third lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalrg  second lowest parts of coefficients of monomials and input, 
 *            space for forward, backward, and cross products;
 *   datalpk  lowest parts of coefficients of monomials and input series, 
 *            space for forward, backward, and cross products;
 *   in1idx   indices of the first input of the addition jobs;
 *   in2idx   indices of the second input of the addition jobs;
 *   outidx   indices of the output of the addition jobs;
 *   dim      the number of coefficients in each series
 *            equals the number of threads in each block.
 *
 * ON RETURN :
 *   datartb  updated highest forward, backward, and cross products;
 *   datarix  updated second highest forward, backward, and cross products;
 *   datarmi  updated third highest forward, backward, and cross products;
 *   datarrg  updated fourth highest forward, backward, and cross products;
 *   datarpk  updated fifth highest forward, backward, and cross products;
 *   dataltb  updated fifth lowest forward, backward, and cross products;
 *   datalix  updated fourth lowest forward, backward, and cross products;
 *   datalmi  updated third lowest forward, backward, and cross products;
 *   datalrg  updated second lowest forward, backward, and cross products;
 *   datalpk  updated lowest forward, backward, and cross products. */

void convoluted_data10_to_output
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions have been computed, but no additions.
 *   This function is only for testing purposes.
 *
 * ON ENTRY :
 *   datartb    highest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   datarix    second highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datarmi    third highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datarrg    fourth highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datarpk    fifth highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   dataltb    fifth lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datalix    fourth lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datalmi    third lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datalrg    second lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products;
 *   datalpk    lowest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products;
 *   outputrtb  has space allocated for dim+1 series of degree deg;
 *   outputrix  has space allocated for dim+1 series of degree deg;
 *   outputrmi  has space allocated for dim+1 series of degree deg;
 *   outputrrg  has space allocated for dim+1 series of degree deg;
 *   outputrpk  has space allocated for dim+1 series of degree deg;
 *   outputltb  has space allocated for dim+1 series of degree deg;
 *   outputlix  has space allocated for dim+1 series of degree deg;
 *   outputlmi  has space allocated for dim+1 series of degree deg;
 *   outputlrg  has space allocated for dim+1 series of degree deg;
 *   outputlpk  has space allocated for dim+1 series of degree deg;
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
 *   outputrtb  has the highest parts of derivatives and the value,
 *              outputrtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrtb[dim] contains the value of the polynomial;
 *   outputrix  has the second highest parts of derivatives and the value,
 *              outputrix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrix[dim] contains the value of the polynomial;
 *   outputrmi  has the third highest parts of derivatives and the value,
 *              outputrmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrmi[dim] contains the value of the polynomial;
 *   outputrrg  has the fourth highest parts of derivatives and the value,
 *              outputrrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrrg[dim] contains the value of the polynomial;
 *   outputrpk  has the fifth highest parts of derivatives and the value,
 *              outputrpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrpk[dim] contains the value of the polynomial;
 *   outputltb  has the fifth lowest parts of derivatives and the value,
 *              outputltb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputltb[dim] contains the value of the polynomial;
 *   outputlix  has the fourth lowest parts of derivatives and the value,
 *              outputlix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlix[dim] contains the value of the polynomial;
 *   outputlmi  has the third lowest parts of derivatives and the value,
 *              outputlmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlmi[dim] contains the value of the polynomial;
 *   outputlrg  has the second lowest parts of derivatives and the value,
 *              outputlrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlrg[dim] contains the value of the polynomial;
 *   outputlpk  has the lowest parts of derivatives and the value,
 *              outputlpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlpk[dim] contains the value of the polynomial. */

void added_data10_to_output
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose=true );
/*
 * DESCRIPTION :
 *   Extracts the data computed on the device to the output.
 *   All convolutions and all additions have been computed.
 *
 * ON ENTRY :
 *   datartb    highest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datarix    second highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datarmi    third highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datarrg    fourth highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datarpk    fifth highest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   dataltb    fifth lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datalix    fourth lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datalmi    third lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datalrg    second lowest parts of coefficients of monomials and input, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   datalpk    lowest parts of coefficients of monomials and input series, 
 *              space for forward, backward, and cross products,
 *              with the accumulated additions;
 *   outputrtb  has space allocated for dim+1 series of degree deg;
 *   outputrix  has space allocated for dim+1 series of degree deg;
 *   outputrmi  has space allocated for dim+1 series of degree deg;
 *   outputrrg  has space allocated for dim+1 series of degree deg;
 *   outputrpk  has space allocated for dim+1 series of degree deg;
 *   outputltb  has space allocated for dim+1 series of degree deg;
 *   outputlix  has space allocated for dim+1 series of degree deg;
 *   outputlmi  has space allocated for dim+1 series of degree deg;
 *   outputlrg  has space allocated for dim+1 series of degree deg;
 *   outputlpk  has space allocated for dim+1 series of degree deg;
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
 *   outputrtb  has the highest parts of derivatives and the value,
 *              outputrtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrtb[dim] contains the value of the polynomial;
 *   outputrix  has the second highest parts of derivatives and the value,
 *              outputrix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrix[dim] contains the value of the polynomial;
 *   outputrmi  has the third highest parts of derivatives and the value,
 *              outputrmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrmi[dim] contains the value of the polynomial;
 *   outputrrg  has the fourth highest parts of derivatives and the value,
 *              outputrrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrrg[dim] contains the value of the polynomial;
 *   outputrpk  has the fifth highest parts of derivatives and the value,
 *              outputrpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrpk[dim] contains the value of the polynomial;
 *   outputltb  has the fifth lowest parts of derivatives and the value,
 *              outputltb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputltb[dim] contains the value of the polynomial;
 *   outputlix  has the fourth lowest parts of derivatives and the value,
 *              outputlix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlix[dim] contains the value of the polynomial;
 *   outputlmi  has the third lowest parts of derivatives and the value,
 *              outputlmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlmi[dim] contains the value of the polynomial;
 *   outputlrg  has the second lowest parts of derivatives and the value,
 *              outputlrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlrg[dim] contains the value of the polynomial;
 *   outputlpk  has the lowest parts of derivatives and the value,
 *              outputlpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlpk[dim] contains the value of the polynomial. */

void GPU_dbl10_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
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
 *   BS          number of threads in a block, must equal deg + 1;
 *   dim         total number of variables;
 *   nbr         number of monomials, excluding the constant term;
 *   deg         truncation degree of the series;
 *   nvr         nvr[k] holds the number of variables in monomial k;
 *   idx         idx[k] has as many indices as the value of nvr[k],
 *               idx[k][i] defines the place of the i-th variable,
 *               with input values in input[idx[k][i]];
 *   cstrtb      highest parts of constant coefficient series;
 *   cstrix      second highest parts of constant coefficient series;
 *   cstrmi      third highest parts of constant coefficient series;
 *   cstrrg      fourth highest parts of constant coefficient series;
 *   cstrpk      fifth highest parts of constant coefficient series;
 *   cstltb      fifth lowest parts of constant coefficient series;
 *   cstlix      fourth lowest parts of constant coefficient series;
 *   cstlmi      third lowest parts of constant coefficient series;
 *   cstlrg      second lowest parts of constant coefficient series;
 *   cstlpk      lowest parts of constant coefficient series;
 *   cffrtb      cffrtb[k] has deg+1 doubles for the highest parts
 *               of the coefficient series of monomial k;
 *   cffrix      cffrix[k] has deg+1 doubles for the second highest parts
 *               of the coefficient series of monomial k;
 *   cffrmi      cffrmi[k] has deg+1 doubles for the third highest parts
 *               of the coefficient series of monomial k;
 *   cffrrg      cffrrg[k] has deg+1 doubles for the fourth highest parts
 *               of the coefficient series of monomial k;
 *   cffrpk      cffrpk[k] has deg+1 doubles for the fifth highest parts
 *               of the coefficient series of monomial k;
 *   cffltb      cffltb[k] has deg+1 doubles for the fifth lowest parts
 *               of the coefficient series of monomial k;
 *   cfflix      cfflix[k] has deg+1 doubles for the fourth lowest parts
 *               of the coefficient series of monomial k;
 *   cfflmi      cfflmi[k] has deg+1 doubles for the third lowest parts
 *               of the coefficient series of monomial k;
 *   cfflrg      cfflrg[k] has deg+1 doubles for the second lowest parts
 *               of the coefficient series of monomial k;
 *   cfflpk      cfflpk[k] has deg+1 doubles for the lowest parts
 *               of the coefficient series of monomial k;
 *   inputrtb    has the highest parts of the power series
 *               for all variables in the polynomial;
 *   inputrix    has the second highest parts of the power series
 *               for all variables in the polynomial;
 *   inputrmi    has the third highest parts of the power series
 *               for all variables in the polynomial;
 *   inputrrg    has the fourth highest parts of the power series
 *               for all variables in the polynomial;
 *   inputrpk    has the fifth highest parts of the power series
 *               for all variables in the polynomial;
 *   inputltb    has the fifth lowest parts of the power series
 *               for all variables in the polynomial;
 *   inputlix    has the fourth lowest parts of the power series
 *               for all variables in the polynomial;
 *   inputlmi    has the third lowest parts of the power series
 *               for all variables in the polynomial;
 *   inputlrg    has the second lowest parts of the power series
 *               for all variables in the polynomial;
 *   inputlpk    has the lowest parts of the power series
 *               for all variables in the polynomial;
 *   outputrtb   has space allocated for dim+1 series of degree deg;
 *   outputrix   has space allocated for dim+1 series of degree deg;
 *   outputrmi   has space allocated for dim+1 series of degree deg;
 *   outputrrg   has space allocated for dim+1 series of degree deg;
 *   outputrpk   has space allocated for dim+1 series of degree deg;
 *   outputltb   has space allocated for dim+1 series of degree deg;
 *   outputlix   has space allocated for dim+1 series of degree deg;
 *   outputlmi   has space allocated for dim+1 series of degree deg;
 *   outputlrg   has space allocated for dim+1 series of degree deg;
 *   outputlpk   has space allocated for dim+1 series of degree deg;
 *   cnvjobs     convolution jobs organized in layers;
 *   addjobs     addition jobs organized in layers;
 *   verbose     if true, then extra output about the setup is written.
 *
 * ON RETURN :
 *   outputrtb   has the highest parts of derivatives and the value,
 *               outputrtb[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputrtb[dim] contains the value of the polynomial;
 *   outputrix   has the second highest parts of derivatives and the value,
 *               outputrix[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputrix[dim] contains the value of the polynomial;
 *   outputrmi   has the third highest parts of derivatives and the value,
 *               outputrmi[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputrmi[dim] contains the value of the polynomial;
 *   outputrrg   has the fourth highest parts of derivatives and the value,
 *               outputrrg[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputrrg[dim] contains the value of the polynomial;
 *   outputrpk   has the fifth highest parts of derivatives and the value,
 *               outputrpk[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputrpk[dim] contains the value of the polynomial;
 *   outputltb   has the fifth lowest parts of derivatives and the value,
 *               outputltb[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputltb[dim] contains the value of the polynomial;
 *   outputlix   has the fourth lowest parts of derivatives and the value,
 *               outputlix[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputlix[dim] contains the value of the polynomial;
 *   outputlmi   has the third lowest parts of derivatives and the value,
 *               outputlmi[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputlmi[dim] contains the value of the polynomial;
 *   outputlrg   has the second lowest parts of derivatives and the value,
 *               outputlrg[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputlrg[dim] contains the value of the polynomial;
 *   outputlpk   has the lowest parts of derivatives and the value,
 *               outputlpk[k], for k from 0 to dim-1, contains the
 *               derivative with respect to the variable k;
 *               outputlpk[dim] contains the value of the polynomial;
 *   cnvlapms    is the elapsed time spent by all convolution kernels,
 *               expressed in milliseconds;
 *   addlapms    is the elapsed time spent by all addition kernels,
 *               expressed in milliseconds;
 *   elapsedms   is the elapsed time spent by all kernels,
 *               expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *               (excluding memory copies) in seconds. */

#endif
