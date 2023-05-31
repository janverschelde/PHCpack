/* The file dbl4_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in quad double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl4_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl4_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl4_polynomials_host_h__
#define __dbl4_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl4_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in quad double precision.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
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
 *   forwardhihi  is work space for the highest doubles of nvr forward 
 *                products, forwardhihi[k] can store deg+1 doubles;
 *   forwardlohi  is work space for the second highest doubles of nvr
 *                forward  products, forwardlohi[k] can store deg+1 doubles;
 *   forwardhilo  is work space for the second lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   forwardlolo  is work space for the lowest doubles of nvr
 *                forward products, forwardlolo[k] can store deg+1 doubles;
 *   backwardhihi is work space for the highest doubles of nvr-1 backward
 *                products, backwardhihi[k] can store deg+1 doubles;
 *   backwardlohi is work space for the second highest doubles of nvr-1
 *                backward products, backwardlohi[k] can store deg+1 doubles;
 *   backwardhilo is work space for the second lowest doubles of nvr-1
 *                backward products, backwardhilo[k] can store deg+1 doubles;
 *   backwardlolo is work space for the lowest doubles of nvr- backward
 *                products, backwardlolo[k] can store deg+1 doubles;
 *   crosshihi    is work space for the highest doubles of nvr-2 cross
 *                products, crosshihi[k] can store deg+1 doubles;
 *   crosslohi    is work space for the second highest doubles of nvr-2
 *                cross products, crosslohi[k] can store deg+1 doubles;
 *   crosshilo    is work space for the second lowest doubles of nvr-2
 *                cross products, crosshilo[k] can store for deg+1 doubles.
 *   crosslolo    is work space for the lowest doubles of nvr-2 cross
 *                products, crosslolo[k] can store for deg+1 doubles.
 *   verbose      if true, writes one line to screen for every convolution.
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
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial. */

void CPU_cmplx4_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihi, double **backwardrelohi,
   double **backwardrehilo, double **backwardrelolo,
   double **backwardimhihi, double **backwardimlohi,
   double **backwardimhilo, double **backwardimlolo,
   double **crossrehihi, double **crossrelohi,
   double **crossrehilo, double **crossrelolo,
   double **crossimhihi, double **crossimlohi,
   double **crossimhilo, double **crossimlolo, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for complex coefficients in quad double precision.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
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
 *   forwardrehihi  has space for the highest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrehihi[k] has space for deg+1 doubles;
 *   forwardrelohi  has space for the second highest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrelohi[k] has space for deg+1 doubles;
 *   forwardrehilo  has space for the second lowest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrehilo[k] has space for deg+1 doubles;
 *   forwardrelolo  has space for the lowest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrelolo[k] has space for deg+1 doubles;
 *   forwardimhihi  has space for the highest doubles of the imaginary parts
 *                  of all nvr forward products,
 *                  forwardimhihi[k] has space for deg+1 doubles;
 *   forwardimlohi  has space for the second highest doubles of 
 *                  the imaginary parts of all nvr forward products,
 *                  forwardimlohi[k] has space for deg+1 doubles;
 *   forwardimhilo  has space for the second lowest doubles of the
 *                  imaginary parts of all nvr forward products,
 *                  forwardimhilo[k] has space for deg+1 doubles;
 *   forwardimlolo  has space for the lowest doubles of the imaginary parts
 *                  of all nvr forward products,
 *                  forwardimlolo[k] has space for deg+1 doubles;
 *   backwardrehihi has space for the highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrehihi[k] has space for deg+1 doubles;
 *   backwardrelohi has space for the second highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrelohi[k] has space for deg+1 doubles;
 *   backwardrehilo has space for the second lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrehilo[k] has space for deg+1 doubles;
 *   backwardrelolo has space for the lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrelolo[k] has space for deg+1 doubles;
 *   backwardimhihi has space for the highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *                  backwardimhihi[k] has space for deg+1 doubles;
 *   backwardimlohi has space for the second highest doubles
 *                  of the imaginary parts of all nvr-2 backward products;
 *                  backwardimlohi[k] has space for deg+1 doubles;
 *   backwardimhilo has space for the second lowest doubles
 *                  of the imaginary parts of all nvr-2 backward products;
 *                  backwardimhilo[k] has space for deg+1 doubles;
 *   backwardimlolo has space for the lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *                  backwardimlolo[k] has space for deg+1 doubles;
 *   crossrehihi    has space for the highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossrehihi[k] has space for deg+1 doubles;
 *   crossrelohi    has space for the second highest doubles
 *                  of the real parts of all nvr-2 cross products;
 *                  crossrelohi[k] has space for deg+1 doubles;
 *   crossrehilo    has space for the second lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossimhilo[k] has space for deg+1 doubles;
 *   crossrelolo    has space for the lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossimlolo[k] has space for deg+1 doubles;
 *   crossimhihi    has space for the highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *                  crossimhihi[k] has space for deg+1 doubles;
 *   crossimlohi    has space for the second highest doubles
 *                  of the imaginary parts of all nvr-2 cross products;
 *                  crossimlohi[k] has space for deg+1 doubles;
 *   crossimhilo    has space for the second lowest doubles 
 *                  of the imaginary parts of all nvr-2 cross products;
 *                  crossimhilo[k] has space for deg+1 doubles;
 *   crossimlolo    has space for the lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *                  crossimlolo[k] has space for deg+1 doubles;
 *   verbose        if true, writes one line to screen for every convolution.
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
 *                  output[dim] contains the value of the polynomial. */

void CPU_dbl4_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
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
 *   vrblvl       is the verbose level.
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
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

void CPU_cmplx4_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *   Evaluates and differentiates the polynomial.
 *
 * ON ENTRY :
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
 *   vrblvl         is the verbose level.
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
 *                  output[dim] contains the value of the polynomial. */

void CPU_dbl4_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   nvr          number of variables in the monomial;
 *   idx          indices to the variables in the monomial;
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
 *   forwardhihi  is work space for the highest doubles of nvr forward 
 *                products, forwardhihi[k] can store deg+1 doubles;
 *   forwardlohi  is work space for the second highest doubles of nvr
 *                forward  products, forwardlohi[k] can store deg+1 doubles;
 *   forwardhilo  is work space for the second lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   forwardlolo  is work space for the lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   backwardhihi is work space for the highest doubles of nvr-1 backward
 *                products, backwardhihi[k] can store deg+1 doubles;
 *   backwardlohi is work space for the second highest doubles of nvr-1
 *                backward products, backwardlohi[k] can store deg+1 doubles;
 *   backwardhilo is work space for the second lowest doubles of nvr-1
 *                backward products, backwardhilo[k] can store deg+1 doubles;
 *   backwardlolo is work space for the lowest doubles of nvr-1 backward
 *                products, backwardlolo[k] can store deg+1 doubles;
 *   crosshihi    is work space for the highest doubles of nvr-2 cross
 *                products, crosshihi[k] can store deg+1 doubles;
 *   crosslohi    is work space for the second highest doubles of nvr-2
 *                cross products, crosslohi[k] can store deg+1 doubles;
 *   crosshilo    is work space for the second lowest doubles of nvr-2
 *                cross products, crosshilo[k] can store for deg+1 doubles;
 *   crosslolo    is work space for the lowest doubles of nvr-2 cross
 *                products, crosslolo[k] can store for deg+1 doubles;
 *   job          defines the convolution job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihi  are the updated highest parts of forward products;
 *   forwardlohi  are the updated second highest parts of forward products;
 *   forwardhilo  are the updated second lowest parts of forward products;
 *   forwardlolo  are the updated lowest parts forward products;
 *   backwardhihi are the updated highest parts of backward products;
 *   backwardlohi are the updated second highest parts of backward products;
 *   backwardhilo are the updated second lowest parts of backward products;
 *   backwardlolo are the updated lowest parts backward products;
 *   crosshihi    are the updated highest parts of cross products;
 *   crosslohi    are the updated second highest parts of cross products;
 *   crosshilo    are the updated second lowest parts of cross products;
 *   crosslolo    are the updated lowest parts cross products. */

void CPU_cmplx4_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehihi, double *cffrelohi,
   double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi,
   double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **forwardrehihi, double **forwardrelohi,
   double **forwardrehilo, double **forwardrelolo,
   double **forwardimhihi, double **forwardimlohi,
   double **forwardimhilo, double **forwardimlolo,
   double **backwardrehihi, double **backwardrelohi,
   double **backwardrehilo, double **backwardrelolo,
   double **backwardimhihi, double **backwardimlohi,
   double **backwardimhilo, double **backwardimlolo,
   double **crossrehihi, double **crossrelohi,
   double **crossrehilo, double **crossrelolo,
   double **crossimhihi, double **crossimlohi,
   double **crossimhilo, double **crossimlolo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   nvr            number of variables in the monomial;
 *   idx            indices to the variables in the monomial;
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
 *   inputrehilo    has the second lowest doubles of the real part
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
 *   forwardrehihi  has space for the highest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrehihi[k] has space for deg+1 doubles;
 *   forwardrelohi  has space for the second highest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrelohi[k] has space for deg+1 doubles;
 *   forwardrehilo  has space for the second lowest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrehilo[k] has space for deg+1 doubles;
 *   forwardrelolo  has space for the lowest doubles of the real parts
 *                  of all nvr forward products,
 *                  forwardrelolo[k] has space for deg+1 doubles;
 *   forwardimhihi  has space for the highest doubles of the imaginary parts
 *                  of all nvr forward products,
 *                  forwardimhihi[k] has space for deg+1 doubles;
 *   forwardimlohi  has space for the second highest doubles of 
 *                  the imaginary parts of all nvr forward products,
 *                  forwardimlohi[k] has space for deg+1 doubles;
 *   forwardimhilo  has space for the second lowest doubles of the
 *                  imaginary parts of all nvr forward products,
 *                  forwardimhilo[k] has space for deg+1 doubles;
 *   forwardimlolo  has space for the lowest doubles of the imaginary parts
 *                  of all nvr forward products,
 *                  forwardimlolo[k] has space for deg+1 doubles;
 *   backwardrehihi has space for the highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrehihi[k] has space for deg+1 doubles;
 *   backwardrelohi has space for the second highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrelohi[k] has space for deg+1 doubles;
 *   backwardrehilo has space for the second lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrehilo[k] has space for deg+1 doubles;
 *   backwardrelolo has space for the lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *                  backwardrelolo[k] has space for deg+1 doubles;
 *   backwardimhihi has space for the highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *                  backwardimhihi[k] has space for deg+1 doubles;
 *   backwardimlohi has space for the second highest doubles
 *                  of the imaginary parts of all nvr-2 backward products;
 *                  backwardimlohi[k] has space for deg+1 doubles;
 *   backwardimhilo has space for the second lowest doubles
 *                  of the imaginary parts of all nvr-2 backward products;
 *                  backwardimhilo[k] has space for deg+1 doubles;
 *   backwardimlolo has space for the lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *                  backwardimlolo[k] has space for deg+1 doubles;
 *   crossrehihi    has space for the highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossrehihi[k] has space for deg+1 doubles;
 *   crossrelohi    has space for the second highest doubles
 *                  of the real parts of all nvr-2 cross products;
 *                  crossrelohi[k] has space for deg+1 doubles;
 *   crossrehilo    has space for the second lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossimhilo[k] has space for deg+1 doubles;
 *   crossrelolo    has space for the lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *                  crossimlolo[k] has space for deg+1 doubles;
 *   crossimhihi    has space for the highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *                  crossimhihi[k] has space for deg+1 doubles;
 *   crossimlohi    has space for the second highest doubles
 *                  of the imaginary parts of all nvr-2 cross products;
 *                  crossimlohi[k] has space for deg+1 doubles;
 *   crossimhilo    has space for the second lowest doubles 
 *                  of the imaginary parts of all nvr-2 cross products;
 *                  crossimhilo[k] has space for deg+1 doubles;
 *   crossimlolo    has space for the lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *                  crossimlolo[k] has space for deg+1 doubles;
 *   job            defines the convolution job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrehihi  are the updated highest doubles of the real parts
 *                  of the forward products;
 *   forwardrelohi  are the updated second highest doubles of the real parts
 *                  of the forward products;
 *   forwardrehilo  are the updated second lowest doubles of the real parts
 *                  of the forward products;
 *   forwardrelolo  are the updated lowest doubles of the real parts
 *                  of the forward products;
 *   forwardimhihi  are the updated highest doubles of the imaginary parts
 *                  of the forward products;
 *   forwardimlohi  are the updated second highest doubles of the
 *                  imaginary parts of the forward products;
 *   forwardimhilo  are the updated second lowest doubles of the imaginary
 *                  parts of the forward products;
 *   forwardimlolo  are the updated lowest doubles of the imaginary parts
 *                  of the forward products;
 *   backwardrehihi are the updated highest doubles of the real parts 
 *                  of the backward products;
 *   backwardrelohi are the updated second highest doubles of the real parts 
 *                  of the backward products;
 *   backwardrehilo are the updated second lowest doubles of the real parts 
 *                  of the backward products;
 *   backwardrelolo are the updated lowest doubles of the real parts 
 *                  of the backward products;
 *   backwardimhihi are the updated highest doubles of the imaginary parts 
 *                  of the backward products;
 *   backwardimlohi are the updated second highest doubles of
 *                  the imaginary parts of the backward products;
 *   backwardimhilo are the updated second lowest doubles of
 *                  the imaginary parts of the backward products;
 *   backwardimlolo are the updated lowest doubles of the imaginary parts 
 *                  of the backward products;
 *   crossrehihi    are the updated highest doubles of the real parts
 *                  of the cross products;
 *   crossrelohi    are the updated second highest doubles of the real parts
 *                  of the cross products;
 *   crossrehilo    are the updated second lowest doubles of the real parts
 *                  of the cross products;
 *   crossrelolo    are the updated lowest doubles of the real parts
 *                  of the cross products;
 *   crossimhihi    are the updated highest doubles of the imaginary parts
 *                  of the cross products;
 *   crossimlohi    are the updated second highest doubles of
 *                  the imaginary parts of the cross products;
 *   crossimhilo    are the updated second lowest doubles of
 *                  the imaginary parts of the cross products;
 *   crossimlolo    are the updated lowest doubles of the imaginary parts
 *                  of the cross products. */

void CPU_dbl4_add_job
 ( int deg,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job, on real data.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
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
 *   forwardhihi  are all highest parts of computed forward products,
 *   forwardlohi  are all second highest parts of computed forward products,
 *   forwardhilo  are all second lowest parts of computed forward products,
 *   forwardlolo  are all lowest parts of computed forward products,
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products;
 *   job          defines the addition job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihi  are the updated highest parts of forward products;
 *   forwardlohi  are the updated second highest parts of forward products;
 *   forwardhilo  are the updated second lowest parts of forward products;
 *   forwardlolo  are the updated lowest parts forward products;
 *   backwardhihi are the updated highest parts of backward products;
 *   backwardlohi are the updated second highest parts of backward products;
 *   backwardhilo are the updated second lowest parts of backward products;
 *   backwardlolo are the updated lowest parts backward products;
 *   crosshihi    are the updated highest parts of cross products;
 *   crosslohi    are the updated second highest parts of cross products;
 *   crosshilo    are the updated second lowest parts of cross products;
 *   crosslolo    are the updated lowest parts cross products. */

void CPU_cmplx4_add_job
 ( int deg, double *cstrehihi, double *cstrelohi,
            double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo, 
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo, 
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
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
 *   forwardrehihi  computed highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelohi  computed second highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrehilo  computed second lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelolo  computed lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardimhihi  computed highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlohi  computed second highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimhilo  computed second lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlolo  computed lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   backwardrehihi computed highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelohi computed second highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrehilo computed second lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelolo computed lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardimhihi computed highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlohi computed second highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimhilo computed second lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlolo computed lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   crossrehihi    computed highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelohi    computed second highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrehilo    computed second lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelolo    computed lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossimhihi    computed highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlohi    computed second highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimhilo    computed second lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlolo    computed lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   job            defines the addition job;
 *   verbose        if true, then is verbose.

 * ON RETURN :
 *   forwardrehihi  are the updated highest doubles of the real parts
 *                  of the forward products;
 *   forwardrelohi  are the updated second highest doubles of the real parts
 *                  of the forward products;
 *   forwardrehilo  are the updated second lowest doubles of the real parts
 *                  of the forward products;
 *   forwardrelolo  are the updated lowest doubles of the real parts
 *                  of the forward products;
 *   forwardimhihi  are the updated highest doubles of the imaginary parts
 *                  of the forward products;
 *   forwardimlohi  are the updated second highest doubles of the
 *                  imaginary parts of the forward products;
 *   forwardimhilo  are the updated second lowest doubles of the
 *                  imaginary parts of the forward products;
 *   forwardimlolo  are the updated lowest doubles of the imaginary parts
 *                  of the forward products;
 *   backwardrehihi are the updated highest doubles of the real parts 
 *                  of the backward products;
 *   backwardrelohi are the updated second highest doubles of the real parts 
 *                  of the backward products;
 *   backwardrehilo are the updated second lowest doubles of the real parts 
 *                  of the backward products;
 *   backwardrelolo are the updated lowest doubles of the real parts 
 *                  of the backward products;
 *   backwardimhihi are the updated highest doubles of the imaginary parts 
 *                  of the backward products;
 *   backwardimlohi are the updated second highest doubles of the
 *                  imaginary parts of the backward products;
 *   backwardimhilo are the updated second lowest doubles of the
 *                  imaginary parts of the backward products;
 *   backwardimlolo are the updated lowest doubles of the imaginary parts 
 *                  of the backward products;
 *   crossrehihi    are the updated highest doubles of the real parts
 *                  of the cross products;
 *   crossrelohi    are the updated second highest doubles of the real parts
 *                  of the cross products;
 *   crossrehilo    are the updated second lowest doubles of the real parts
 *                  of the cross products;
 *   crossrelolo    are the updated lowest doubles of the real parts
 *                  of the cross products;
 *   crossimhihi    are the updated highest doubles of the
 *                  imaginary parts of the cross products;
 *   crossimlohi    are the updated second highest doubles of the
 *                  imaginary parts of the cross products;
 *   crossimhilo    are the updated lowest doubles of the imaginary parts
 *                  of the cross products;
 *   crossimlolo    are the updated lowest doubles of the imaginary parts
 *                  of the cross products. */

void CPU_dbl4_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
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
 *   forwardhihi  are all highest parts of computed forward products;
 *   forwardlohi  are all second highest parts of computed forward products;
 *   forwardhilo  are all second lowest parts of computed forward products;
 *   forwardlolo  are all lowest parts of computed forward products;
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products.
 *
 * ON RETURN :
 *   outputhihi   highest parts of the values and all derivatives;
 *   outputlohi   second highest parts of the values and all derivatives;
 *   outputhilo   second lowest parts of the values and all derivatives;
 *   outputlolo   lowest parts of the values and all derivatives. */

void CPU_cmplx4_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo,
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo,
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
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
 *   forwardrehihi  computed highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelohi  computed second highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrehilo  computed second lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelolo  computed lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardimhihi  computed highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlohi  computed second highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimhilo  computed second lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlolo  computed lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   backwardrehihi computed high doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelohi computed second highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrehilo computed second lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelolo computed lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardimhihi computed highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlohi computed second highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimhilo computed second lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlolo computed lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   crossrehihi    computed high doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelohi    computed second highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrehilo    computed second lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelolo    computed lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossimhihi    computed highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlohi    computed second highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimhilo    computed second lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlolo    computed lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products.
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
 *                  output[dim] contains the value of the polynomial. */

void CPU_dbl4_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
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
 *   forwardhihi  are all highest parts of computed forward products;
 *   forwardlohi  are all second highest parts of computed forward products;
 *   forwardhilo  are all second lowest parts of computed forward products;
 *   forwardlolo  are all lowest parts of computed forward products;
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products;
 *   jobs         defines the addition jobs;
 *   verbose      if true, then output is written.
 *
 * ON RETURN :
 *   outputhihi   highest parts of the values and all derivatives;
 *   outputlohi   second highest parts of the values and all derivatives;
 *   outputhilo   second lowest parts of the values and all derivatives;
 *   outputlolo   lowest parts of the values and all derivatives. */

void CPU_cmplx4_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrehihi, double ***forwardrelohi,
   double ***forwardrehilo, double ***forwardrelolo,
   double ***forwardimhihi, double ***forwardimlohi,
   double ***forwardimhilo, double ***forwardimlolo,
   double ***backwardrehihi, double ***backwardrelohi,
   double ***backwardrehilo, double ***backwardrelolo,
   double ***backwardimhihi, double ***backwardimlohi,
   double ***backwardimhilo, double ***backwardimlolo,
   double ***crossrehihi, double ***crossrelohi,
   double ***crossrehilo, double ***crossrelolo,
   double ***crossimhihi, double ***crossimlohi,
   double ***crossimhilo, double ***crossimlolo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs to the final output,
 *   on real data.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
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
 *   forwardrehihi  computed highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelohi  computed second highest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrehilo  computed second lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardrelolo  computed lowest doubles of the real parts
 *                  of all nvr forward products;
 *   forwardimhihi  computed highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlohi  computed second highest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimhilo  computed second lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   forwardimlolo  computed lowest doubles of the imaginary parts
 *                  of all nvr forward products;
 *   backwardrehihi computed high doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelohi computed second highest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrehilo computed second lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardrelolo computed lowest doubles of the real parts 
 *                  of all nvr-2 backward products;
 *   backwardimhihi computed highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlohi computed second highest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimhilo computed second lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   backwardimlolo computed lowest doubles of the imaginary parts 
 *                  of all nvr-2 backward products;
 *   crossrehihi    computed high doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelohi    computed second highest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrehilo    computed second lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossrelolo    computed lowest doubles of the real parts
 *                  of all nvr-2 cross products;
 *   crossimhihi    computed highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlohi    computed second highest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimhilo    computed second lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products;
 *   crossimlolo    computed lowest doubles of the imaginary parts
 *                  of all nvr-2 cross products.
 *   jobs           defines the addition jobs;
 *   verbose        if true, then output is written.
 *
 * ON RETURN :
 *   outputrehihi   are the highest doubles of the real parts
 *                  of the values and all derivatives;
 *   outputrelohi   are the second highest doubles of the real parts
 *                  of the values and all derivatives;
 *   outputrehilo   are the second lowest doubles of the real parts
 *                  of the values and all derivatives;
 *   outputrelolo   are the lowest doubles of the real parts
 *                  of the values and all derivatives;
 *   outputimhihi   are the highest doubles of the imaginary parts
 *                  of the values and all derivatives;
 *   outputimlohi   are the second highest doubles of the imaginary parts
 *                  of the values and all derivatives;
 *   outputimhilo   are the second lowest doubles of the imaginary parts
 *                  of the values and all derivatives;
 *   outputimlolo   are the lowest doubles of the imaginary parts
 *                  of the values and all derivatives. */

void CPU_dbl4_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
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
 *   vrblvl       is the verbose level.
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
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

void CPU_cmplx4_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs to the final output,
 *   on real data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
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
 *   vrblvl         is the verbose level.
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
 *   elapsedsec     is the elapsed time in seconds. */

#endif
