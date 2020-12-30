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

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
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
 *   verbose  if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

void CPU_dbl_conv_job
 ( int deg, int nvr, int *idx, double *cff, double **input,
   double **forward, double **backward, double **cross,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
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

void CPU_dbl_add_job
 ( int deg, double *cst, double **cff,
   double ***forward, double ***backward, double ***cross,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job.
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

void CPU_dbl_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
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

void CPU_dbl_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   double ***forward, double ***backward, double ***cross,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
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

void CPU_dbl_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose=false );
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
 *   verbose  if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

#endif
