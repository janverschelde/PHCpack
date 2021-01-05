/* The file dbl10_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in penta double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl10_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl10_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl10_polynomials_host_h__
#define __dbl10_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl10_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in deca double precision.
 *
 * ON ENTRY :
 *   dim         total number of variables;
 *   nbr         number of monomials, excluding the constant term;
 *   deg         truncation degree of the series;
 *   nvr         nvr[k] holds the number of variables in monomial k;
 *   idx         idx[k] has as many indices as the value of nvr[k],
 *               idx[k][i] defines the place of the i-th variable,
 *               with input values in input[idx[k][i]];
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
 *   forwardrtb  is work space for the highest doubles of nvr forward 
 *               products, forwardrtb[k] can store deg+1 doubles;
 *   forwardrix  is work space for the second highest doubles of nvr
 *               forward  products, forwardrix[k] can store deg+1 doubles;
 *   forwardrmi  is work space for the third highest doubles of nvr
 *               forward  products, forwardrmi[k] can store deg+1 doubles;
 *   forwardrrg  is work space for the fourth highest doubles of nvr
 *               forward  products, forwardrrg[k] can store deg+1 doubles;
 *   forwardrpk  is work space for the fifth highest doubles of nvr
 *               forward  products, forwardrpk[k] can store deg+1 doubles;
 *   forwardltb  is work space for the fifth lowest doubles of nvr
 *               forward  products, forwardltb[k] can store deg+1 doubles;
 *   forwardlix  is work space for the fourth lowest doubles of nvr
 *               forward  products, forwardlix[k] can store deg+1 doubles;
 *   forwardlmi  is work space for the third lowest doubles of nvr
 *               forward  products, forwardlmi[k] can store deg+1 doubles;
 *   forwardlrg  is work space for the second lowest doubles of nvr
 *               forward products, forwardlrg[k] can store deg+1 doubles;
 *   forwardlpk  is work space for the lowest doubles of nvr
 *               forward products, forwardlpk[k] can store deg+1 doubles;
 *   backwardrtb is work space for the highest doubles of nvr-1 backward
 *               products, backwardtb[k] can store deg+1 doubles;
 *   backwardrix is work space for the second highest doubles of nvr-1
 *               backward products, backwardrix[k] can store deg+1 doubles;
 *   backwardrmi is work space for the third highest doubles of nvr-1
 *               backward products, backwardrmi[k] can store deg+1 doubles;
 *   backwardrrg is work space for the fourth highest doubles of nvr-1
 *               backward products, backwardrrg[k] can store deg+1 doubles;
 *   backwardrpk is work space for the fifth highest doubles of nvr-1
 *               backward products, backwardrpk[k] can store deg+1 doubles;
 *   backwardltb is work space for the fifth lowest doubles of nvr-1
 *               backward products, backwardltb[k] can store deg+1 doubles;
 *   backwardlix is work space for the fourth lowest doubles of nvr-1
 *               backward products, backwardlix[k] can store deg+1 doubles;
 *   backwardlmi is work space for the third lowest doubles of nvr-1
 *               backward products, backwardlmi[k] can store deg+1 doubles;
 *   backwardlrg is work space for the second lowest doubles of nvr-1
 *               backward products, backwardlrg[k] can store deg+1 doubles;
 *   backwardlpk is work space for the lowest doubles of nvr-1 backward
 *               products, backwardlpk[k] can store deg+1 doubles;
 *   crossrtb    is work space for the highest doubles of nvr-2 cross
 *               products, crossrtb[k] can store deg+1 doubles;
 *   crossrix    is work space for the second highest doubles of nvr-2
 *               cross products, crossrix[k] can store deg+1 doubles;
 *   crossrmi    is work space for the third highest doubles of nvr-2
 *               cross products, crossrmi[k] can store deg+1 doubles;
 *   crossrrg    is work space for the fourth highest doubles of nvr-2
 *               cross products, crossrrg[k] can store deg+1 doubles;
 *   crossrpk    is work space for the fifth highest doubles of nvr-2
 *               cross products, crossrpk[k] can store deg+1 doubles;
 *   crossltb    is work space for the fifth lowest doubles of nvr-2
 *               cross products, crossltb[k] can store deg+1 doubles;
 *   crosslix    is work space for the fourth lowest doubles of nvr-2
 *               cross products, crosslix[k] can store deg+1 doubles;
 *   crosslmi    is work space for the third lowest doubles of nvr-2
 *               cross products, crosslmi[k] can store deg+1 doubles;
 *   crosslrg    is work space for the second lowest doubles of nvr-2
 *               cross products, crosslrg[k] can store for deg+1 doubles;
 *   crosslpk    is work space for the lowest doubles of nvr-2 cross
 *               products, crosslpk[k] can store for deg+1 doubles;
 *   verbose     if true, writes one line to screen for every convolution.
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
 *               outputlpk[dim] contains the value of the polynomial. */

void CPU_dbl10_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
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
 *   verbose     if true, writes one line to screen for every convolution.
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
 *   elapsedsec  is the elapsed time in seconds. */

void CPU_dbl10_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
 *
 * ON ENTRY :
 *   deg         degree of the series;
 *   nvr         number of variables in the monomial;
 *   idx         indices to the variables in the monomial;
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
 *   forwardrtb  is work space for the highest doubles of nvr forward 
 *               products, forwardrtb[k] can store deg+1 doubles;
 *   forwardrix  is work space for the second highest doubles of nvr
 *               forward  products, forwardrix[k] can store deg+1 doubles;
 *   forwardrmi  is work space for the third highest doubles of nvr
 *               forward  products, forwardrmi[k] can store deg+1 doubles;
 *   forwardrrg  is work space for the fourth highest doubles of nvr
 *               forward  products, forwardrrg[k] can store deg+1 doubles;
 *   forwardrpk  is work space for the fifth highest doubles of nvr
 *               forward  products, forwardrpk[k] can store deg+1 doubles;
 *   forwardltb  is work space for the fifth lowest doubles of nvr
 *               forward  products, forwardltb[k] can store deg+1 doubles;
 *   forwardlix  is work space for the fourth lowest doubles of nvr
 *               forward  products, forwardlix[k] can store deg+1 doubles;
 *   forwardlmi  is work space for the third lowest doubles of nvr
 *               forward  products, forwardlmi[k] can store deg+1 doubles;
 *   forwardlrg  is work space for the second lowest doubles of nvr
 *               forward products, forwardlrg[k] can store deg+1 doubles;
 *   forwardlpk  is work space for the lowest doubles of nvr
 *               forward products, forwardlpk[k] can store deg+1 doubles;
 *   backwardrtb is work space for the highest doubles of nvr-1 backward
 *               products, backwardtb[k] can store deg+1 doubles;
 *   backwardrix is work space for the second highest doubles of nvr-1
 *               backward products, backwardrix[k] can store deg+1 doubles;
 *   backwardrmi is work space for the third highest doubles of nvr-1
 *               backward products, backwardrmi[k] can store deg+1 doubles;
 *   backwardrrg is work space for the fourth highest doubles of nvr-1
 *               backward products, backwardrrg[k] can store deg+1 doubles;
 *   backwardrpk is work space for the fifth highest doubles of nvr-1
 *               backward products, backwardrpk[k] can store deg+1 doubles;
 *   backwardltb is work space for the fifth lowest doubles of nvr-1
 *               backward products, backwardltb[k] can store deg+1 doubles;
 *   backwardlix is work space for the fourth lowest doubles of nvr-1
 *               backward products, backwardlix[k] can store deg+1 doubles;
 *   backwardlmi is work space for the third lowest doubles of nvr-1
 *               backward products, backwardlmi[k] can store deg+1 doubles;
 *   backwardlrg is work space for the second lowest doubles of nvr-1
 *               backward products, backwardlrg[k] can store deg+1 doubles;
 *   backwardlpk is work space for the lowest doubles of nvr-1 backward
 *               products, backwardlpk[k] can store deg+1 doubles;
 *   crossrtb    is work space for the highest doubles of nvr-2 cross
 *               products, crossrtb[k] can store deg+1 doubles;
 *   crossrix    is work space for the second highest doubles of nvr-2
 *               cross products, crossrix[k] can store deg+1 doubles;
 *   crossrmi    is work space for the third highest doubles of nvr-2
 *               cross products, crossrmi[k] can store deg+1 doubles;
 *   crossrrg    is work space for the fourth highest doubles of nvr-2
 *               cross products, crossrrg[k] can store deg+1 doubles;
 *   crossrpk    is work space for the fifth highest doubles of nvr-2
 *               cross products, crossrpk[k] can store deg+1 doubles;
 *   crossltb    is work space for the fifth lowest doubles of nvr-2
 *               cross products, crossltb[k] can store deg+1 doubles;
 *   crosslix    is work space for the fourth lowest doubles of nvr-2
 *               cross products, crosslix[k] can store deg+1 doubles;
 *   crosslmi    is work space for the third lowest doubles of nvr-2
 *               cross products, crosslmi[k] can store deg+1 doubles;
 *   crosslrg    is work space for the second lowest doubles of nvr-2
 *               cross products, crosslrg[k] can store for deg+1 doubles;
 *   crosslpk    is work space for the lowest doubles of nvr-2 cross
 *               products, crosslpk[k] can store for deg+1 doubles;
 *   job         defines the convolution job;
 *   verbose     if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrtb  are the updated highest parts of forward products;
 *   forwardrix  are the updated second highest parts of forward products;
 *   forwardrmi  are the updated third highest parts of forward products;
 *   forwardrrg  are the updated fourth highest parts of forward products;
 *   forwardrpk  are the updated fifth highest parts of forward products;
 *   forwardltb  are the updated fifth lowest parts of forward products;
 *   forwardlix  are the updated fourth lowest parts of forward products;
 *   forwardlmi  are the updated third lowest parts of forward products;
 *   forwardlrg  are the updated second lowest parts of forward products;
 *   forwardlpk  are the updated lowest parts forward products;
 *   backwardrtb are the updated highest parts of backward products;
 *   backwardrix are the updated second highest parts of backward products;
 *   backwardrmi are the updated third highest parts of backward products;
 *   backwardrrg are the updated fourth highest parts of backward products;
 *   backwardrpk are the updated fifth highest parts of backward products;
 *   backwardltb are the updated fifth lowest parts of backward products;
 *   backwardlmi are the updated fourth lowest parts of backward products;
 *   backwardlmi are the updated third lowest parts of backward products;
 *   backwardlrg are the updated second lowest parts of backward products;
 *   backwardlpk are the updated lowest parts backward products;
 *   crossrtb    are the updated highest parts of cross products;
 *   crossrix    are the updated second highest parts of cross products;
 *   crossrmi    are the updated third highest parts of cross products;
 *   crossrrg    are the updated fourth highest parts of cross products;
 *   crossrpk    are the updated fifth highest parts of cross products;
 *   crossltb    are the updated fifth lowest parts of cross products;
 *   crosslix    are the updated fourth lowest parts of cross products;
 *   crosslmi    are the updated third lowest parts of cross products;
 *   crosslrg    are the updated second lowest parts of cross products;
 *   crosslpk    are the updated lowest parts cross products. */

void CPU_dbl10_add_job
 ( int deg,
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job.
 *
 * ON ENTRY :
 *   deg         degree of the series;
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
 *   forwardrtb  are all highest parts of computed forward products;
 *   forwardrix  are all second highest parts of computed forward products;
 *   forwardrmi  are all third highest parts of computed forward products;
 *   forwardrrg  are all fourth highest parts of computed forward products;
 *   forwardrpk  are all fifth highest parts of computed forward products;
 *   forwardltb  are all fifth lowest parts of computed forward products;
 *   forwardlix  are all fourth lowest parts of computed forward products;
 *   forwardlmi  are all third lowest parts of computed forward products;
 *   forwardlrg  are all second lowest parts of computed forward products;
 *   forwardlpk  are all lowest parts of computed forward products;
 *   backwardrtb are all highest parts of computed backward products;
 *   backwardrix are all second highest parts of computed backward products;
 *   backwardrmi are all third highest parts of computed backward products;
 *   backwardrrg are all fourth highest parts of computed backward products;
 *   backwardrpk are all fifth highest parts of computed backward products;
 *   backwardltb are all fifth lowest parts of computed backward products;
 *   backwardlix are all fourth lowest parts of computed backward products;
 *   backwardlmi are all third lowest parts of computed backward products;
 *   backwardlrg are all second lowest parts of computed backward products;
 *   backwardlpk are all lowest parts of computed backward products;
 *   crossrtb    are all highest parts of computed cross products;
 *   crossrix    are all second highest parts of computed cross products;
 *   crossrmi    are all third highest parts of computed cross products;
 *   crossrrg    are all fourth highest parts of computed cross products;
 *   crossrpk    are all fifth highest parts of computed cross products;
 *   crossltb    are all fifth lowest parts of computed cross products;
 *   crosslix    are all fourth lowest parts of computed cross products;
 *   crosslmi    are all third lowest parts of computed cross products;
 *   crosslrg    are all second lowest parts of computed cross products;
 *   crosslpk    are all lowest parts of computed cross products;
 *   job         defines the addition job;
 *   verbose     if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrtb  are the updated highest parts of forward products;
 *   forwardrix  are the updated second highest parts of forward products;
 *   forwardrmi  are the updated third highest parts of forward products;
 *   forwardrrg  are the updated fourth highest parts of forward products;
 *   forwardrpk  are the updated fifth highest parts of forward products;
 *   forwardltb  are the updated fifth lowest parts of forward products;
 *   forwardlix  are the updated fourth lowest parts of forward products;
 *   forwardlmi  are the updated third lowest parts of forward products;
 *   forwardlrg  are the updated second lowest parts of forward products;
 *   forwardlpk  are the updated lowest parts forward products;
 *   backwardrtb are the updated highest parts of backward products;
 *   backwardrix are the updated second highest parts of backward products;
 *   backwardrmi are the updated third highest parts of backward products;
 *   backwardrrg are the updated fourth highest parts of backward products;
 *   backwardrpk are the updated fifth highest parts of backward products;
 *   backwardltb are the updated fifth lowest parts of backward products;
 *   backwardlmi are the updated fourth lowest parts of backward products;
 *   backwardlmi are the updated third lowest parts of backward products;
 *   backwardlrg are the updated second lowest parts of backward products;
 *   backwardlpk are the updated lowest parts backward products;
 *   crossrtb    are the updated highest parts of cross products;
 *   crossrix    are the updated second highest parts of cross products;
 *   crossrmi    are the updated third highest parts of cross products;
 *   crossrrg    are the updated fourth highest parts of cross products;
 *   crossrpk    are the updated fifth highest parts of cross products;
 *   crossltb    are the updated fifth lowest parts of cross products;
 *   crosslix    are the updated fourth lowest parts of cross products;
 *   crosslmi    are the updated third lowest parts of cross products;
 *   crosslrg    are the updated second lowest parts of cross products;
 *   crosslpk    are the updated lowest parts cross products. */

void CPU_dbl10_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
 *
 * ON ENTRY :
 *   dim         total number of variables;
 *   nbr         number of monomials, excluding the constant term;
 *   deg         degree of the series;
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
 *   forwardrtb  are all highest parts of computed forward products;
 *   forwardrix  are all second highest parts of computed forward products;
 *   forwardrmi  are all third highest parts of computed forward products;
 *   forwardrrg  are all fourth highest parts of computed forward products;
 *   forwardrpk  are all fifth highest parts of computed forward products;
 *   forwardltb  are all fifth lowest parts of computed forward products;
 *   forwardlix  are all fourth lowest parts of computed forward products;
 *   forwardlmi  are all third lowest parts of computed forward products;
 *   forwardlrg  are all second lowest parts of computed forward products;
 *   forwardlpk  are all lowest parts of computed forward products;
 *   backwardrtb are all highest parts of computed backward products;
 *   backwardrix are all second highest parts of computed backward products;
 *   backwardrmi are all third highest parts of computed backward products;
 *   backwardrrg are all fourth highest parts of computed backward products;
 *   backwardrpk are all fifth highest parts of computed backward products;
 *   backwardltb are all fifth lowest parts of computed backward products;
 *   backwardlix are all fourth lowest parts of computed backward products;
 *   backwardlmi are all third lowest parts of computed backward products;
 *   backwardlrg are all second lowest parts of computed backward products;
 *   backwardlpk are all lowest parts of computed backward products;
 *   crossrtb    are all highest parts of computed cross products;
 *   crossrix    are all second highest parts of computed cross products;
 *   crossrmi    are all third highest parts of computed cross products;
 *   crossrrg    are all fourth highest parts of computed cross products;
 *   crossrpk    are all fifth highest parts of computed cross products;
 *   crossltb    are all fifth lowest parts of computed cross products;
 *   crosslix    are all fourth lowest parts of computed cross products;
 *   crosslmi    are all third lowest parts of computed cross products;
 *   crosslrg    are all second lowest parts of computed cross products;
 *   crosslpk    are all lowest parts of computed cross products.
 *
 * ON RETURN :
 *   outputrtb   highest parts of the values and all derivatives;
 *   outputrix   second highest parts of the values and all derivatives;
 *   outputrmi   third highest parts of the values and all derivatives;
 *   outputrrg   fourth highest parts of the values and all derivatives;
 *   outputrpk   fifth highest parts of the values and all derivatives;
 *   outputltb   fifth lowest parts of the values and all derivatives;
 *   outputlix   fourth lowest parts of the values and all derivatives;
 *   outputlmi   third lowest parts of the values and all derivatives;
 *   outputlrg   second lowest parts of the values and all derivatives;
 *   outputlpk   lowest parts of the values and all derivatives. */

void CPU_dbl10_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrtb, double ***forwardrix, double ***forwardrmi,
   double ***forwardrrg, double ***forwardrpk,
   double ***forwardltb, double ***forwardlix, double ***forwardlmi,
   double ***forwardlrg, double ***forwardlpk,
   double ***backwardrtb, double ***backwardrix, double ***backwardrmi,
   double ***backwardrrg, double ***backwardrpk, 
   double ***backwardltb, double ***backwardlix, double ***backwardlmi,
   double ***backwardlrg, double ***backwardlpk, 
   double ***crossrtb, double ***crossrix, double ***crossrmi,
   double ***crossrrg, double ***crossrpk,
   double ***crossltb, double ***crosslix, double ***crosslmi,
   double ***crosslrg, double ***crosslpk,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
 *
 * ON ENTRY :
 *   dim         total number of variables;
 *   nbr         number of monomials, excluding the constant term;
 *   deg         degree of the series;
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
 *   forwardrtb  are all highest parts of computed forward products;
 *   forwardrix  are all second highest parts of computed forward products;
 *   forwardrmi  are all third highest parts of computed forward products;
 *   forwardrrg  are all fourth highest parts of computed forward products;
 *   forwardrpk  are all fifth highest parts of computed forward products;
 *   forwardltb  are all fifth lowest parts of computed forward products;
 *   forwardlix  are all fourth lowest parts of computed forward products;
 *   forwardlmi  are all third lowest parts of computed forward products;
 *   forwardlrg  are all second lowest parts of computed forward products;
 *   forwardlpk  are all lowest parts of computed forward products;
 *   backwardrtb are all highest parts of computed backward products;
 *   backwardrix are all second highest parts of computed backward products;
 *   backwardrmi are all third highest parts of computed backward products;
 *   backwardrrg are all fourth highest parts of computed backward products;
 *   backwardrpk are all fifth highest parts of computed backward products;
 *   backwardltb are all fifth lowest parts of computed backward products;
 *   backwardlix are all fourth lowest parts of computed backward products;
 *   backwardlmi are all third lowest parts of computed backward products;
 *   backwardlrg are all second lowest parts of computed backward products;
 *   backwardlpk are all lowest parts of computed backward products;
 *   crossrtb    are all highest parts of computed cross products;
 *   crossrix    are all second highest parts of computed cross products;
 *   crossrmi    are all third highest parts of computed cross products;
 *   crossrrg    are all fourth highest parts of computed cross products;
 *   crossrpk    are all fifth highest parts of computed cross products;
 *   crossltb    are all fifth lowest parts of computed cross products;
 *   crosslix    are all fourth lowest parts of computed cross products;
 *   crosslmi    are all third lowest parts of computed cross products;
 *   crosslrg    are all second lowest parts of computed cross products;
 *   crosslpk    are all lowest parts of computed cross products;
 *   jobs        defines the addition jobs;
 *   verbose     if true, then output is written.
 *
 * ON RETURN :
 *   outputrtb   highest parts of the values and all derivatives;
 *   outputrix   second highest parts of the values and all derivatives;
 *   outputrmi   third highest parts of the values and all derivatives;
 *   outputrrg   fourth highest parts of the values and all derivatives;
 *   outputrpk   fifth highest parts of the values and all derivatives;
 *   outputltb   fifth lowest parts of the values and all derivatives;
 *   outputlix   fourth lowest parts of the values and all derivatives;
 *   outputlmi   third lowest parts of the values and all derivatives;
 *   outputlrg   second lowest parts of the values and all derivatives;
 *   outputlpk   lowest parts of the values and all derivatives. */

void CPU_dbl10_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
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
 *   verbose     if true, writes one line to screen for every convolution.
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
 *   elapsedsec  is the elapsed time in seconds. */

#endif
