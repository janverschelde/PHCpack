/* The file dbl_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl_polynomials_host_h__
#define __dbl_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cff, double **input, double **output,
   double **forward, double **backward, double **cross, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cff      cff[k] has deg+1 doubles for the coefficient series
 *            of monomial k;
 *   input    contains the coefficients of the power series
 *            for all variables in the polynomial;
 *   output   space allocated for the value and all derivatives;
 *   forward  contains work space for all nvr forward products,
 *            forward[k] contains space for deg+1 doubles;
 *   backward contains work space for all nvr-2 backward products;
 *            backward[k] contains space for deg+1 doubles;
 *   cross    contains work space for all nvr-2 cross products;
 *            cross[k] contains space for deg+1 doubles;
 *   verbose  if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

void CPU_cmplx_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffre, double **cffim, double **inputre, double **inputim,
   double **outputre, double **outputim,
   double **forwardre, double **forwardim,
   double **backwardre, double **backwardim,
   double **crossre, double **crossim, bool verbose );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for complex coefficients in double precision.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   cffre      cffre[k] has deg+1 doubles for the real parts 
 *              of the coefficient series of monomial k;
 *   cffim      cffim[k] has deg+1 doubles for the imaginary parts
 *              of the coefficient series of monomial k;
 *   inputre    contains the real parts of power series coefficients;
 *   inputim    contains the imaginary parts of power series coefficients;
 *   outputre   has space allocated for the value and all derivatives;
 *   outputim   has space allocated for the value and all derivatives;
 *   forwardre  has work space for all nvr forward products,
 *              forwardre[k] contains space for deg+1 doubles;
 *   forwardim  has work space for all nvr forward products,
 *              forwardim[k] contains space for deg+1 doubles;
 *   backwardre has work space for all nvr-2 backward products;
 *              backwardre[k] contains space for deg+1 doubles;
 *   backwardim has work space for all nvr-2 backward products;
 *              backwardim[k] contains space for deg+1 doubles;
 *   crossre    has work space for all nvr-2 cross products;
 *              crossre[k] contains space for deg+1 doubles;
 *   crossim    has work space for all nvr-2 cross products;
 *              crossim[k] contains space for deg+1 doubles;
 *   verbose    if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputre   contains the real parts of the derivatives and
 *              the value of the polynomial,
 *              outputre[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputre[dim] contains the value of the polynomial;
 *   outputim   contains the imaginary parts of the derivatives and
 *              the of the polynomial,
 *              outputim[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputim[dim] contains the value of the polynomial. */

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *   Evaluates and differentiates the polynomial.
 *
 * ON ENTRY :
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
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial;
 *   elapsedsec is the elapsed time in seconds. */

void CPU_cmplx_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim, double **cffre, double **cffim,
   double **inputre, double **inputim, double **outputre, double **outputim,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *   Evaluates and differentiates the polynomial.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cstre    deg+1 doubles for the real parts of the constant;
 *   cstim    deg+1 doubles for the imaginary parts of the constant;
 *   cffre    cffre[k] has deg+1 doubles for the real parts 
 *            of the coefficient series of monomial k;
 *   cffim    cffim[k] has deg+1 doubles for the imaginary parts 
 *            of the coefficient series of monomial k;
 *   inputre  has the real parts of the power series coefficients;
 *   inputim  has the imaginary parts of the power series coefficients;
 *   outputre has space allocated for the value and all derivatives;
 *   outputim has space allocated for the value and all derivatives;
 *   output   space allocated for the value and all derivatives;
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
 *   elapsedsec is the elapsed time in seconds. */

void CPU_dbl_conv_job
 ( int deg, int nvr, int *idx, double *cff, double **input,
   double **forward, double **backward, double **cross,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job on real data.
 *
 * ON ENTRY :
 *   deg      degree of the series;
 *   nvr      number of variables in the monomial;
 *   idx      indices to the variables in the monomial;
 *   cff      coefficient series of the monomial;
 *   input    power series for all variables in the polynomial;
 *   forward  contains work space for all nvr forward products,
 *            forward[k] contains space for deg+1 doubles;
 *   backward contains work space for all nvr-2 backward products;
 *            backward[k] contains space for deg+1 doubles;
 *   cross    contains work space for all nvr-2 cross products;
 *            cross[k] contains space for deg+1 doubles;
 *   job      defines the convolution job;
 *   verbose  if true, then is verbose.
 *
 * ON RETURN :
 *   forward  are the updated forward products;
 *   backward are the updated backward products;
 *   cross    are the updated cross products. */

void CPU_cmplx_conv_job
 ( int deg, int nvr, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim,
   double **forwardre, double **forwardim,
   double **backwardre, double **backwardim,
   double **crossre, double **crossim,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job on complex data.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   nvr        number of variables in the monomial;
 *   idx        indices to the variables in the monomial;
 *   cffre      real parts of the coefficient series of the monomial;
 *   cffim      imaginary parts of the coefficient series of the monomial;
 *   inputre    real parts of the power series for all variables;
 *   inputim    imaginary parts of the power series for all variables;
 *   forwardre  has work space for all nvr forward products,
 *              forwardre[k] has space for deg+1 doubles;
 *   forwardim  has work space for all nvr forward products,
 *              forwardim[k] has space for deg+1 doubles;
 *   backwardre has work space for all nvr-2 backward products;
 *              backwardre[k] has space for deg+1 doubles;
 *   backwardim has work space for all nvr-2 backward products;
 *              backwardim[k] has space for deg+1 doubles;
 *   crossre    has work space for all nvr-2 cross products;
 *              crossre[k] has space for deg+1 doubles;
 *   crossim    has work space for all nvr-2 cross products;
 *              crossim[k] has space for deg+1 doubles;
 *   job        defines the convolution job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardre  are the updated real parts of the forward products;
 *   forwardim  are the updated imaginary parts of the forward products;
 *   backwardre are the updated real parts of the backward products;
 *   backwardim are the updated imaginary parts of the backward products;
 *   crossre    are the updated real parts of the cross products;
 *   crossim    are the updated imaginary parts of the cross products. */

void CPU_dbl_add_job
 ( int deg, double *cst, double **cff,
   double ***forward, double ***backward, double ***cross,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job on real data.
 *
 * ON ENTRY :
 *   deg      degree of the series;
 *   cst      constant coefficient series of the polynmial;
 *   cff      coefficients of the monomials;
 *   forward  all computed forward products,
 *   backward all computed backward products;
 *   cross    all computed cross products;
 *   job      defines the addition job;
 *   verbose  if true, then is verbose.
 *
 * ON RETURN :
 *   forward  are the updated forward products;
 *   backward are the updated backward products;
 *   cross    are the updated cross products. */

void CPU_cmplx_add_job
 ( int deg, double *cstre, double *cstim,
   double **cffre, double **cffim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job on complex data.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   cstre      real parts of the constant coefficient series;
 *   cstim      imaginary parts of the constant coefficient series;
 *   cffre      real parts of the coefficients of the monomials;
 *   cffim      imaginary parts of the coefficients of the monomials;
 *   forwardre  has all computed real parts of the forward products,
 *   forwardim  has all computed imaginary parts of the forward products,
 *   backwardre has all computed real parts of the backward products;
 *   backwardim has all computed imaginary parts of the backward products;
 *   crossre    has all computed real parts of the cross products;
 *   crossim    has all computed imaginary parts of the cross products;
 *   job        defines the addition job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardre  are the updated real parts of the forward products;
 *   forwardim  are the updated imaginary parts of the forward products;
 *   backwardre are the updated real parts of the backward products;
 *   backwardim are the updated imaginary parts of the backward products;
 *   crossre    are the updated real parts of the cross products;
 *   crossim    are the updated imaginary parts of the cross products. */

void CPU_dbl_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output,
 *   for real data.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cst      constant coefficient series of the polynmial;
 *   cff      coefficients of the monomials;
 *   forward  are all computed forward products,
 *   backward are all computed backward products;
 *   cross    are all computed cross products.
 *
 * ON RETURN :
 *   output   the values and all derivatives. */

void CPU_cmplx_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output,
 *   for complex data.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   cstre      real parts of the constant coefficient series;
 *   cstim      imaginary parts of the constant coefficient series;
 *   cffre      real parts of the coefficients of the monomials;
 *   cffim      imaginary parts of the coefficients of the monomials;
 *   forwardre  are all computed real parts of the forward products,
 *   forwardim  are all computed imaginary parts of the forward products,
 *   backwardre are all computed real parts of the backward products;
 *   backwardim are all computed imaginary parts of the backward products;
 *   crossre    are all computed real parts of the cross products.
 *   crossim    are all computed imaginary parts of the cross products.
 *
 * ON RETURN :
 *   outputre   the real parts of the values and all derivatives;
 *   outputim   the imaginary parts of the values and all derivatives. */

void CPU_dbl_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs on real data.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cst      constant coefficient series of the polynmial;
 *   cff      coefficients of the monomials;
 *   forward  are all computed forward products,
 *   backward are all computed backward products;
 *   cross    are all computed cross products;
 *   jobs     defines the addition jobs;
 *   verbose  if true, then output is written.
 *
 * ON RETURN :
 *   output   the values and all derivatives. */

void CPU_cmplx_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   double ***forwardre, double ***forwardim,
   double ***backwardre, double ***backwardim,
   double ***crossre, double ***crossim,
   AdditionJobs jobs, bool verbose );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs on complex data.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   cstre      real parts of the constant coefficient series;
 *   cstim      imaginary parts of the constant coefficient series;
 *   cffre      real parts of the coefficients of the monomials;
 *   cffim      imaginary parts of the coefficients of the monomials;
 *   forwardre  are all computed real parts of the forward products,
 *   forwardim  are all computed imaginary parts of the forward products,
 *   backwardre are all computed real parts of the backward products;
 *   backwardim are all computed imaginary parts of the backward products;
 *   crossre    are all computed real parts of the cross products;
 *   crossim    are all computed imaginary parts of the cross products;
 *   jobs       defines the addition jobs;
 *   verbose    if true, then output is written.
 *
 * ON RETURN :
 *   outputre   the real parts of the values and all derivatives;
 *   outputim   the imaginary parts of the values and all derivatives. */

void CPU_dbl_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
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
 *   elapsedsec is the elapsed time in seconds. */

void CPU_cmplx_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstre, double *cstim,
   double **cffre, double **cffim,
   double **inputre, double **inputim,
   double **outputre, double **outputim,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cstre    deg+1 doubles for the real parts of the constant;
 *   cstim    deg+1 doubles for the imaginary parts of the constant;
 *   cffre    cffre[k] has deg+1 doubles for the real parts
 *            of the coefficient series of monomial k;
 *   cffim    cffim[k] has deg+1 doubles for the imaginary part
 *            of the coefficient series of monomial k;
 *   inputre  has the real parts of the coefficients of the power series
 *            for all variables in the polynomial;
 *   inputim  has the imaginary parts of the coefficients of the power series
 *            for all variables in the polynomial;
 *   outputre has space allocated for the value and all derivatives;
 *   outputim has space allocated for the value and all derivatives;
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
 *   elapsedsec is the elapsed time in seconds. */

#endif
