/* The file dbl8_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in octo double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl8_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl8_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl8_polynomials_host_h__
#define __dbl8_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl8_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double **forwardhihihi, double **forwardhilohi,
   double **forwardhihilo, double **forwardhilolo,
   double **forwardlohihi, double **forwardlolohi,
   double **forwardlohilo, double **forwardlololo,
   double **backwardhihihi, double **backwardhilohi,
   double **backwardhihilo, double **backwardhilolo,
   double **backwardlohihi, double **backwardlolohi,
   double **backwardlohilo, double **backwardlololo,
   double **crosshihihi, double **crosshilohi,
   double **crosshihilo, double **crosshilolo,
   double **crosslohihi, double **crosslolohi,
   double **crosslohilo, double **crosslololo,
   bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in octo double precision.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   forwardhihihi  is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardlohihi  is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardhilohi  is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardlolohi  is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardhihilo  is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardlohilo  is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardhilolo  is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardlololo  is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   backwardhihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   crosshihihi    is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosslohihi    is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosshilohi    is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosslolohi    is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosshihilo    is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosslohilo    is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosshilolo    is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosslololo    is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   verbose        if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi   has the second highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo   has the third highest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo   has the fourth highest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi   has the fourth lowest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi   has the third lowest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo   has the second lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlololo[dim] contains the value of the polynomial. */

void CPU_dbl8_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   csthihihi      highest parts of constant coefficient series;
 *   csthilohi      second higest parts of constant coefficient series;
 *   csthihilo      third higest parts of constant coefficient series;
 *   csthilolo      fourth higest parts of constant coefficient series;
 *   cstlohihi      fourth lowest parts of constant coefficient series;
 *   cstlolohi      third lowest parts of constant coefficient series;
 *   cstlohilo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg;
 *   verbose        if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi   has the second highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo   has the third highest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo   has the fourth highest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi   has the fourth lowest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi   has the third lowest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo   has the second lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlololo[dim] contains the value of the polynomial. */

void CPU_dbl8_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihihi, double *cffhilohi,
   double *cffhihilo, double *cffhilolo,
   double *cfflohihi, double *cfflolohi,
   double *cfflohilo, double *cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo,
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo,
   double **forwardhihihi, double **forwardhilohi,
   double **forwardhihilo, double **forwardhilolo,
   double **forwardlohihi, double **forwardlolohi,
   double **forwardlohilo, double **forwardlololo,
   double **backwardhihihi, double **backwardhilohi,
   double **backwardhihilo, double **backwardhilolo,
   double **backwardlohihi, double **backwardlolohi,
   double **backwardlohilo, double **backwardlololo,
   double **crosshihihi, double **crosshilohi,
   double **crosshihilo, double **crosshilolo,
   double **crosslohihi, double **crosslolohi,
   double **crosslohilo, double **crosslololo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   nvr            number of variables in the monomial;
 *   idx            indices to the variables in the monomial;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg;
 *   forwardhihihi  is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardlohihi  is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardhilohi  is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardlolohi  is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardhihilo  is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardlohilo  is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardhilolo  is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardlololo  is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   backwardhihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardlololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   crosshihihi    is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosslohihi    is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosshilohi    is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosslolohi    is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crosshihilo    is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosslohilo    is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosshilolo    is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crosslololo    is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   job            defines the convolution job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihihi  are the updated highest parts of forward products;
 *   forwardhilohi  are the updated second highest parts of forward products;
 *   forwardhihilo  are the updated third highest parts of forward products;
 *   forwardhilolo  are the updated fourth highest parts of forward products;
 *   forwardlohihi  are the updated fourth lowest parts of forward products;
 *   forwardlohilo  are the updated third lowest parts of forward products;
 *   forwardlolohi  are the updated second lowest parts of forward products;
 *   forwardlololo  are the updated lowest parts forward products;
 *   backwardhihihi are the updated highest parts of backward products;
 *   backwardhilohi are the updated second highest parts of backward products;
 *   backwardhihilo are the updated third highest parts of backward products;
 *   backwardhilolo are the updated fourth highest parts of backward products;
 *   backwardlohihi are the updated fourth lowest parts of backward products;
 *   backwardlolohi are the updated third lowest parts of backward products;
 *   backwardlohilo are the updated second lowest parts of backward products;
 *   backwardlololo are the updated lowest parts backward products;
 *   crosshihihi    are the updated highest parts of cross products;
 *   crosshilohi    are the updated second highest parts of cross products;
 *   crosshihilo    are the updated third highest parts of cross products;
 *   crosshilolo    are the updated fourth highest parts of cross products;
 *   crosslohihi    are the updated fourth lowest parts of cross products;
 *   crosslolohi    are the updated third lowest parts of cross products;
 *   crosslohilo    are the updated second lowest parts of cross products;
 *   crosslololo    are the updated lowest parts cross products. */

void CPU_dbl8_add_job
 ( int deg,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   csthihihi      highest parts of constant coefficient series;
 *   csthilohi      second higest parts of constant coefficient series;
 *   csthihilo      third higest parts of constant coefficient series;
 *   csthilolo      fourth higest parts of constant coefficient series;
 *   cstlohihi      fourth lowest parts of constant coefficient series;
 *   cstlolohi      third lowest parts of constant coefficient series;
 *   cstlohilo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   forwardhihihi  are all highest parts of computed forward products;
 *   forwardhilohi  are all second highest parts of computed forward products;
 *   forwardhihilo  are all third highest parts of computed forward products;
 *   forwardhilolo  are all fourth highest parts of computed forward products;
 *   forwardlohihi  are all fourth lowest parts of computed forward products;
 *   forwardlohilo  are all third lowest parts of computed forward products;
 *   forwardlolohi  are all second lowest parts of computed forward products;
 *   forwardlololo  are all lowest parts of computed forward products;
 *   backwardhihihi are all highest parts of computed backward products;
 *   backwardhilohi are all second highest parts of computed backward products;
 *   backwardhihilo are all third highest parts of computed backward products;
 *   backwardhilolo are all fourth highest parts of computed backward products;
 *   backwardlohilo are all second lowest parts of computed backward products;
 *   backwardlololo are all lowest parts of computed backward products;
 *   crosshihihi    are all highest parts of computed cross products;
 *   crosshilohi    are all second highest parts of computed cross products;
 *   crosshihilo    are all third highest parts of computed cross products;
 *   crosshilolo    are all fourth highest parts of computed cross products;
 *   crosslohihi    are all fourth lowest parts of computed cross products;
 *   crosslohilo    are all third lowest parts of computed cross products;
 *   crosslolohi    are all second lowest parts of computed cross products;
 *   crosslololo    are all lowest parts of computed cross products;
 *   job            defines the addition job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihihi  are the updated highest parts of forward products;
 *   forwardhilohi  are the updated second highest parts of forward products;
 *   forwardhihilo  are the updated third highest parts of forward products;
 *   forwardhilolo  are the updated fourth highest parts of forward products;
 *   forwardlohihi  are the updated fourth lowest parts of forward products;
 *   forwardlohilo  are the updated third lowest parts of forward products;
 *   forwardlolohi  are the updated second lowest parts of forward products;
 *   forwardlololo  are the updated lowest parts forward products;
 *   backwardhihihi are the updated highest parts of backward products;
 *   backwardhilohi are the updated second highest parts of backward products;
 *   backwardhihilo are the updated third highest parts of backward products;
 *   backwardhilolo are the updated fourth highest parts of backward products;
 *   backwardlohihi are the updated fourth lowest parts of backward products;
 *   backwardlolohi are the updated third lowest parts of backward products;
 *   backwardlohilo are the updated second lowest parts of backward products;
 *   backwardlololo are the updated lowest parts backward products;
 *   crosshihihi    are the updated highest parts of cross products;
 *   crosshilohi    are the updated second highest parts of cross products;
 *   crosshihilo    are the updated third highest parts of cross products;
 *   crosshilolo    are the updated fourth highest parts of cross products;
 *   crosslohihi    are the updated fourth lowest parts of cross products;
 *   crosslolohi    are the updated third lowest parts of cross products;
 *   crosslohilo    are the updated second lowest parts of cross products;
 *   crosslololo    are the updated lowest parts cross products. */

void CPU_dbl8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   csthihihi      highest parts of constant coefficient series;
 *   csthilohi      second higest parts of constant coefficient series;
 *   csthihilo      third higest parts of constant coefficient series;
 *   csthilolo      fourth higest parts of constant coefficient series;
 *   cstlohihi      fourth lowest parts of constant coefficient series;
 *   cstlolohi      third lowest parts of constant coefficient series;
 *   cstlohilo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   forwardhihihi  are all highest parts of computed forward products;
 *   forwardhilohi  are all second highest parts of computed forward products;
 *   forwardhihilo  are all third highest parts of computed forward products;
 *   forwardhilolo  are all fourth highest parts of computed forward products;
 *   forwardlohihi  are all fourth lowest parts of computed forward products;
 *   forwardlohilo  are all third lowest parts of computed forward products;
 *   forwardlolohi  are all second lowest parts of computed forward products;
 *   forwardlololo  are all lowest parts of computed forward products;
 *   backwardhihihi are all highest parts of computed backward products;
 *   backwardhilohi are all second highest parts of computed backward products;
 *   backwardhihilo are all third highest parts of computed backward products;
 *   backwardhilolo are all fourth highest parts of computed backward products;
 *   backwardlohilo are all second lowest parts of computed backward products;
 *   backwardlololo are all lowest parts of computed backward products;
 *   crosshihihi    are all highest parts of computed cross products;
 *   crosshilohi    are all second highest parts of computed cross products;
 *   crosshihilo    are all third highest parts of computed cross products;
 *   crosshilolo    are all fourth highest parts of computed cross products;
 *   crosslohihi    are all fourth lowest parts of computed cross products;
 *   crosslohilo    are all third lowest parts of computed cross products;
 *   crosslolohi    are all second lowest parts of computed cross products;
 *   crosslololo    are all lowest parts of computed cross products.
 *
 * ON RETURN :
 *   outputhihihi   highest parts of the values and all derivatives;
 *   outputhilohi   second highest parts of the values and all derivatives;
 *   outputhihilo   third highest parts of the values and all derivatives;
 *   outputhilolo   fourth highest parts of the values and all derivatives;
 *   outputlohihi   fourth lowest parts of the values and all derivatives;
 *   outputlolohi   third lowest parts of the values and all derivatives;
 *   outputlohilo   second lowest parts of the values and all derivatives;
 *   outputlololo   lowest parts of the values and all derivatives. */

void CPU_dbl8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardhihihi, double ***forwardhilohi,
   double ***forwardhihilo, double ***forwardhilolo,
   double ***forwardlohihi, double ***forwardlolohi,
   double ***forwardlohilo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardhilohi,
   double ***backwardhihilo, double ***backwardhilolo, 
   double ***backwardlohihi, double ***backwardlolohi,
   double ***backwardlohilo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosshilohi,
   double ***crosshihilo, double ***crosshilolo,
   double ***crosslohihi, double ***crosslolohi,
   double ***crosslohilo, double ***crosslololo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   csthihihi      highest parts of constant coefficient series;
 *   csthilohi      second higest parts of constant coefficient series;
 *   csthihilo      third higest parts of constant coefficient series;
 *   csthilolo      fourth higest parts of constant coefficient series;
 *   cstlohihi      fourth lowest parts of constant coefficient series;
 *   cstlolohi      third lowest parts of constant coefficient series;
 *   cstlohilo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   forwardhihihi  are all highest parts of computed forward products;
 *   forwardhilohi  are all second highest parts of computed forward products;
 *   forwardhihilo  are all third highest parts of computed forward products;
 *   forwardhilolo  are all fourth highest parts of computed forward products;
 *   forwardlohihi  are all fourth lowest parts of computed forward products;
 *   forwardlohilo  are all third lowest parts of computed forward products;
 *   forwardlolohi  are all second lowest parts of computed forward products;
 *   forwardlololo  are all lowest parts of computed forward products;
 *   backwardhihihi are all highest parts of computed backward products;
 *   backwardhilohi are all second highest parts of computed backward products;
 *   backwardhihilo are all third highest parts of computed backward products;
 *   backwardhilolo are all fourth highest parts of computed backward products;
 *   backwardlohilo are all second lowest parts of computed backward products;
 *   backwardlololo are all lowest parts of computed backward products;
 *   jobs           defines the addition jobs;
 *   verbose        if true, then output is written.
 *
 * ON RETURN :
 *   outputhihihi   highest parts of the values and all derivatives;
 *   outputhilohi   second highest parts of the values and all derivatives;
 *   outputhihilo   third highest parts of the values and all derivatives;
 *   outputhilolo   fourth highest parts of the values and all derivatives;
 *   outputlohihi   fourth lowest parts of the values and all derivatives;
 *   outputlolohi   third lowest parts of the values and all derivatives;
 *   outputlohilo   second lowest parts of the values and all derivatives;
 *   outputlololo   lowest parts of the values and all derivatives. */

void CPU_dbl8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   csthihihi      highest parts of constant coefficient series;
 *   csthilohi      second higest parts of constant coefficient series;
 *   csthihilo      third higest parts of constant coefficient series;
 *   csthilolo      fourth higest parts of constant coefficient series;
 *   cstlohihi      fourth lowest parts of constant coefficient series;
 *   cstlolohi      third lowest parts of constant coefficient series;
 *   cstlohilo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg;
 *   cnvjobs        convolution jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   verbose        if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihihi[dim] contains the value of the polynomial;
 *   outputhilohi   has the second highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilohi[dim] contains the value of the polynomial;
 *   outputhihilo   has the third highest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihilo[dim] contains the value of the polynomial;
 *   outputhilolo   has the fourth highest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilolo[dim] contains the value of the polynomial;
 *   outputlohihi   has the fourth lowest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohihi[dim] contains the value of the polynomial;
 *   outputlolohi   has the third lowest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlolohi[dim] contains the value of the polynomial;
 *   outputlohilo   has the second lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlololo[dim] contains the value of the polynomial. */

#endif
