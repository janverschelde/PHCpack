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
   double **forwardhihihi, double **forwardlohihi,
   double **forwardhilohi, double **forwardlolohi,
   double **forwardhihilo, double **forwardlohilo,
   double **forwardhilolo, double **forwardlololo,
   double **backwardhihihi, double **backwardlohihi,
   double **backwardhilohi, double **backwardlolohi,
   double **backwardhihilo, double **backwardlohilo,
   double **backwardhilolo, double **backwardlololo,
   double **crosshihihi, double **crosslohihi,
   double **crosshilohi, double **crosslolohi,
   double **crosshihilo, double **crosslohilo,
   double **crosshilolo, double **crosslololo,
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
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhiloihi[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
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
 *   verbose        if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k,
 *                  the value of the polynomial is at position dim;
 *   outputlohihi   has the second highest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputhilohi   has the third highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputlolohi   has the fourth highest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputhihilo   has the fourth lowest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputlohilo   has the third lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputhilolo   has the second lowest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k.
 *                  the value of the polynomial is at position dim. */

void CPU_cmplx8_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double **forwardrehihihi, double **forwardrelohihi,
   double **forwardrehilohi, double **forwardrelolohi,
   double **forwardrehihilo, double **forwardrelohilo,
   double **forwardrehilolo, double **forwardrelololo,
   double **forwardimhihihi, double **forwardimlohihi,
   double **forwardimhilohi, double **forwardimlolohi,
   double **forwardimhihilo, double **forwardimlohilo,
   double **forwardimhilolo, double **forwardimlololo,
   double **backwardrehihihi, double **backwardrelohihi,
   double **backwardrehilohi, double **backwardrelolohi,
   double **backwardrehihilo, double **backwardrelohilo,
   double **backwardrehilolo, double **backwardrelololo,
   double **backwardimhihihi, double **backwardimlohihi,
   double **backwardimhilohi, double **backwardimlolohi,
   double **backwardimhihilo, double **backwardimlohilo,
   double **backwardimhilolo, double **backwardimlololo,
   double **crossrehihihi, double **crossrelohihi,
   double **crossrehilohi, double **crossrelolohi,
   double **crossrehihilo, double **crossrelohilo,
   double **crossrehilolo, double **crossrelololo,
   double **crossimhihihi, double **crossimlohihi,
   double **crossimhilohi, double **crossimlolohi,
   double **crossimhihilo, double **crossimlohilo,
   double **crossimhilolo, double **crossimlololo,
   bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for complex coefficients in octo double precision.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            truncation degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   cffrehihihi    has the highest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrelohihi    has the second highest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrehilohi    has the third highest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrelolohi    has the fourth highest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrehihilo    has the fourth lowest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrelohilo    has the third lowest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrehilolo    has the second lowest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffrelololo    has the lowest doubles of the real parts
 *                  of the coefficient series of monomial k;
 *   cffimhihihi    has the highest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimlohihi    has the second highest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimhilohi    has the third highest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimlolohi    has the fourth highest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimhihilo    has the fourth lowest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimlohilo    has the third lowest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimhilolo    has the second lowest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   cffimlololo    has the lowest doubles of the imaginary parts
 *                  of the coefficient series of monomial k;
 *   inputrehihihi  has the highest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrelohihi  has the second highest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrehilohi  has the third highest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrelolohi  has the fourth highest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrehihilo  has the fourth lowest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrelohilo  has the third lowest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrehilolo  has the second lowest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputrelololo  has the lowest parts of the real inputs
 *                  for all variables in the polynomial;
 *   inputimhihihi  has the highest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimlohihi  has the second highest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimhilohi  has the third highest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimlolohi  has the fourth highest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimhihilo  has the fourth lowest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimlohilo  has the third lowest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimhilolo  has the second lowest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   inputimlololo  has the lowest parts of the imaginary inputs
 *                  for all variables in the polynomial;
 *   outputrehihihi has space allocated for dim+1 series of degree deg;
 *   outputrelohihi has space allocated for dim+1 series of degree deg;
 *   outputrehilohi has space allocated for dim+1 series of degree deg;
 *   outputrelolohi has space allocated for dim+1 series of degree deg;
 *   outputrehihilo has space allocated for dim+1 series of degree deg;
 *   outputrelohilo has space allocated for dim+1 series of degree deg;
 *   outputrehilolo has space allocated for dim+1 series of degree deg;
 *   outputrelololo has space allocated for dim+1 series of degree deg;
 *   outputimhihihi has space allocated for dim+1 series of degree deg;
 *   outputimlohihi has space allocated for dim+1 series of degree deg;
 *   outputimhilohi has space allocated for dim+1 series of degree deg;
 *   outputimlolohi has space allocated for dim+1 series of degree deg;
 *   outputimhihilo has space allocated for dim+1 series of degree deg;
 *   outputimlohilo has space allocated for dim+1 series of degree deg;
 *   outputimhilolo has space allocated for dim+1 series of degree deg;
 *   outputimlohilo has space allocated for dim+1 series of degree deg;
 *   forwardrehihihi is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrelohihi is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrehilohi is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrelolohi is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrehihilo is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrelohilo is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrehilolo is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardrelololo is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimhihihi is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimlohihi is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimhilohi is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimlolohi is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimhihilo is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimlohilo is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimhilolo is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardimlololo is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   backwardrehihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   crossrehihihi  is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrelohihi    is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrehilohi    is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrelolohi    is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrehihilo    is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrelohilo    is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrehilolo    is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrelololo    is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimhihihi    is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimlohihi    is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimhilohi    is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimlolohi    is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimhilo    is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimhilo    is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimlolo    is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossololo    is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   verbose        if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputrehihihi has the highest parts of derivatives and the value,
 *                  outputrehihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k,
 *                  the value of the polynomial is at position dim;
 *   outputrelohihi has the second highest parts of derivatives and the value,
 *                  outputrelohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrehilohi has the third highest parts of derivatives and the value,
 *                  outputrehilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrelolohi has the fourth highest parts of derivatives and the value,
 *                  outputrelolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrehihilo has the fourth lowest parts of derivatives and the value,
 *                  outputrehihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrelohilo has the third lowest parts of derivatives and the value,
 *                  outputrelohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrehilolo has the second lowest parts of derivatives and the value,
 *                  outputrehilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputrelololo has the lowest parts of derivatives and the value,
 *                  outputrelololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k.
 *                  the value of the polynomial is at position dim;
 *   outputimhihihi has the highest parts of derivatives and the value,
 *                  outputimhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k,
 *                  the value of the polynomial is at position dim;
 *   outputimlohihi has the second highest parts of derivatives and the value,
 *                  outputimlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimhilohi has the third highest parts of derivatives and the value,
 *                  outputimhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimlolohi has the fourth highest parts of derivatives and the value,
 *                  outputimlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimhihilo has the fourth lowest parts of derivatives and the value,
 *                  outputimhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimlohilo has the third lowest parts of derivatives and the value,
 *                  outputimlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimhilolo has the second lowest parts of derivatives and the value,
 *                  outputimhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  the value of the polynomial is at position dim;
 *   outputimlololo has the lowest parts of derivatives and the value,
 *                  outputimlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k.
 *                  the value of the polynomial is at position dim. */

void CPU_dbl8_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, int vrblvl );
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
 *   csthilohi      second highest parts of constant coefficient series;
 *   csthilohi      third highest parts of constant coefficient series;
 *   cstlolohi      fourth highest parts of constant coefficient series;
 *   csthihilo      fourth lowest parts of constant coefficient series;
 *   cstlohilo      third lowest parts of constant coefficient series;
 *   csthilolo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cfflohihi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihihi[dim] contains the value of the polynomial;
 *   outputlohihi   has the second highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilohi[dim] contains the value of the polynomial;
 *   outputhilohi   has the third highest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihilo[dim] contains the value of the polynomial;
 *   outputlolohi   has the fourth highest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilolo[dim] contains the value of the polynomial;
 *   outputhihilo   has the fourth lowest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohihi[dim] contains the value of the polynomial;
 *   outputlohilo   has the third lowest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlolohi[dim] contains the value of the polynomial;
 *   outputhilolo   has the second lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohilo[dim] contains the value of the polynomial;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlololo[dim] contains the value of the polynomial;
 *   elapsedsec     is the elapsed time in seconds. */

void CPU_cmplx8_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, int vrblvl );
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
 *   cstrehihihi    highest parts of constant real coefficients;
 *   cstrehilohi    second highest parts of constant real coefficients;
 *   cstrehilohi    third highest parts of constant real coefficients;
 *   cstrelolohi    fourth highest parts of constant real coefficients;
 *   cstrehihilo    fourth lowest parts of constant real coefficients;
 *   cstrelohilo    third lowest parts of constant real coefficients;
 *   cstrehilolo    second lowest parts of constant real coefficients;
 *   cstrelololo    lowest parts of constant real coefficients;
 *   cstimhihihi    highest parts of constant imag coefficients;
 *   cstimlolohi    second highest parts of constant imag coefficients;
 *   cstimhilohi    third highest parts of constant imag coefficients;
 *   cstimlolohi    fourth highest parts of constant imag coefficients;
 *   cstimhihilo    fourth lowest parts of constant imag coefficients;
 *   cstimlohilo    third lowest parts of constant imag coefficients;
 *   cstimhilolo    second lowest parts of constant imag coefficients;
 *   cstimlololo    lowest parts of constant imag coefficients;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrehilohi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilohi    cffrehihilo[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrehilolo[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrelohihi[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelolohi[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrelohilo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimhilohi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilohi    cffimhihilo[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimhilolo[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimlohihi[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlolohi[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimlohilo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   inputrehihihi  has the highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohihi  has the second highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilohi  has the third highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelolohi  has the fourth highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehihilo  has the fourth lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohilo  has the third lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilolo  has the second lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelololo  has the lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihihi  has the highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohihi  has the second highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilohi  has the third highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlolohi  has the fourth highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihilo  has the fourth lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohilo  has the third lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilolo  has the second lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlololo  has the lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   outputrehihihi has space allocated for dim+1 series of degree deg;
 *   outputrelohihi has space allocated for dim+1 series of degree deg;
 *   outputrehilohi has space allocated for dim+1 series of degree deg;
 *   outputrelolohi has space allocated for dim+1 series of degree deg;
 *   outputrehihilo has space allocated for dim+1 series of degree deg;
 *   outputrelohilo has space allocated for dim+1 series of degree deg;
 *   outputrehilolo has space allocated for dim+1 series of degree deg;
 *   outputrelololo has space allocated for dim+1 series of degree deg;
 *   outputimhihihi has space allocated for dim+1 series of degree deg;
 *   outputimlohihi has space allocated for dim+1 series of degree deg;
 *   outputimhilohi has space allocated for dim+1 series of degree deg;
 *   outputimlolohi has space allocated for dim+1 series of degree deg;
 *   outputimhihilo has space allocated for dim+1 series of degree deg;
 *   outputimlohilo has space allocated for dim+1 series of degree deg;
 *   outputimhilolo has space allocated for dim+1 series of degree deg;
 *   outputimlololo has space allocated for dim+1 series of degree deg;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputrehihihi has the highest parts of derivatives and the value,
 *                  outputrehihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrehihihi[dim] contains the value of the polynomial;
 *   outputrelohihi has the second highest parts of derivatives and the value,
 *                  outputrehilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrehilohi[dim] contains the value of the polynomial;
 *   outputrehilohi has the third highest parts of derivatives and the value,
 *                  outputrehihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrehihilo[dim] contains the value of the polynomial;
 *   outputrelolohi has the fourth highest parts of derivatives and the value,
 *                  outputrehilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrehilolo[dim] contains the value of the polynomial;
 *   outputrehihilo has the fourth lowest parts of derivatives and the value,
 *                  outputrelohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrelohihi[dim] contains the value of the polynomial;
 *   outputrelohilo has the third lowest parts of derivatives and the value,
 *                  outputrelolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrelolohi[dim] contains the value of the polynomial;
 *   outputrehilolo has the second lowest parts of derivatives and the value,
 *                  outputrelohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrelohilo[dim] contains the value of the polynomial;
 *   outputrelololo has the lowest parts of derivatives and the value,
 *                  outputrelololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputrelololo[dim] contains the value of the polynomial;
 *   outputimhihihi has the highest parts of derivatives and the value,
 *                  outputimhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimhihihi[dim] contains the value of the polynomial;
 *   outputimlohihi has the second highest parts of derivatives and the value,
 *                  outputimhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimhilohi[dim] contains the value of the polynomial;
 *   outputimhilohi has the third highest parts of derivatives and the value,
 *                  outputimhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimhihilo[dim] contains the value of the polynomial;
 *   outputimlolohi has the fourth highest parts of derivatives and the value,
 *                  outputimhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimhilolo[dim] contains the value of the polynomial;
 *   outputimhihilo has the fourth lowest parts of derivatives and the value,
 *                  outputimlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimlohihi[dim] contains the value of the polynomial;
 *   outputimlohilo has the third lowest parts of derivatives and the value,
 *                  outputimlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimlolohi[dim] contains the value of the polynomial;
 *   outputimhilolo has the second lowest parts of derivatives and the value,
 *                  outputimlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimlohilo[dim] contains the value of the polynomial;
 *   outputimlololo has the lowest parts of derivatives and the value,
 *                  outputimlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputimlololo[dim] contains the value of the polynomial;
 *   elapsedsec     is the elapsed time in seconds. */

void CPU_dbl8_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **forwardhihihi, double **forwardlohihi,
   double **forwardhilohi, double **forwardlolohi,
   double **forwardhihilo, double **forwardlohilo,
   double **forwardhilolo, double **forwardlololo,
   double **backwardhihihi, double **backwardlohihi,
   double **backwardhilohi, double **backwardlolohi,
   double **backwardhihilo, double **backwardlohilo,
   double **backwardhilolo, double **backwardlololo,
   double **crosshihihi, double **crosslohihi,
   double **crosshilohi, double **crosslolohi,
   double **crosshihilo, double **crosslohilo,
   double **crosshilolo, double **crosslololo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job on real data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   nvr            number of variables in the monomial;
 *   idx            indices to the variables in the monomial;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   inputhihihi    has the highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohihi    has the second highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilohi    has the third highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlolohi    has the fourth highest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhihilo    has the fourth lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlohilo    has the third lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputhilolo    has the second lowest parts of the power series
 *                  for all variables in the polynomial;
 *   inputlololo    has the lowest parts of the power series
 *                  for all variables in the polynomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
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
 *   forwardlohihi  are the updated second highest parts of forward products;
 *   forwardhilohi  are the updated third highest parts of forward products;
 *   forwardlolohi  are the updated fourth highest parts of forward products;
 *   forwardhihilo  are the updated fourth lowest parts of forward products;
 *   forwardlohilo  are the updated third lowest parts of forward products;
 *   forwardhilolo  are the updated second lowest parts of forward products;
 *   forwardlololo  are the updated lowest parts forward products;
 *   backwardhihihi are the updated highest parts of backward products;
 *   backwardlohihi are the updated second highest parts of backward products;
 *   backwardhilohi are the updated third highest parts of backward products;
 *   backwardlolohi are the updated fourth highest parts of backward products;
 *   backwardhihilo are the updated fourth lowest parts of backward products;
 *   backwardlohilo are the updated third lowest parts of backward products;
 *   backwardhilolo are the updated second lowest parts of backward products;
 *   backwardlololo are the updated lowest parts backward products;
 *   crosshihihi    are the updated highest parts of cross products;
 *   crosslohihi    are the updated second highest parts of cross products;
 *   crosshilohi    are the updated third highest parts of cross products;
 *   crosslolohi    are the updated fourth highest parts of cross products;
 *   crosshihilo    are the updated fourth lowest parts of cross products;
 *   crosslohilo    are the updated third lowest parts of cross products;
 *   crosshilolo    are the updated second lowest parts of cross products;
 *   crosslololo    are the updated lowest parts cross products. */

void CPU_cmplx8_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **forwardrehihihi, double **forwardrelohihi,
   double **forwardrehilohi, double **forwardrelolohi,
   double **forwardrehihilo, double **forwardrelohilo,
   double **forwardrehilolo, double **forwardrelololo,
   double **forwardimhihihi, double **forwardimlohihi,
   double **forwardimhilohi, double **forwardimlolohi,
   double **forwardimhihilo, double **forwardimlohilo,
   double **forwardimhilolo, double **forwardimlololo,
   double **backwardrehihihi, double **backwardrelohihi,
   double **backwardrehilohi, double **backwardrelolohi,
   double **backwardrehihilo, double **backwardrelohilo,
   double **backwardrehilolo, double **backwardrelololo,
   double **backwardimhihihi, double **backwardimlohihi,
   double **backwardimhilohi, double **backwardimlolohi,
   double **backwardimhihilo, double **backwardimlohilo,
   double **backwardimhilolo, double **backwardimlololo,
   double **crossrehihihi, double **crossrelohihi,
   double **crossrehilohi, double **crossrelolohi,
   double **crossrehihilo, double **crossrelohilo,
   double **crossrehilolo, double **crossrelololo,
   double **crossimhihihi, double **crossimlohihi,
   double **crossimhilohi, double **crossimlolohi,
   double **crossimhihilo, double **crossimlohilo,
   double **crossimhilolo, double **crossimlololo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job on complex data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   nvr            number of variables in the monomial;
 *   idx            indices to the variables in the monomial;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrelohihi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilohi    cffrehilohi[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrelolohi[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrehihilo[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelohilo[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrehilolo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimlohihi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilohi    cffimhilohi[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimlolohi[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimhihilo[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlohilo[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimhilolo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   inputrehihihi  has the highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohihi  has the second highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilohi  has the third highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelolohi  has the fourth highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehihilo  has the fourth lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohilo  has the third lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilolo  has the second lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelololo  has the lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihihi  has the highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohihi  has the second highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilohi  has the third highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlolohi  has the fourth highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihilo  has the fourth lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohilo  has the third lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilolo  has the second lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlololo  has the lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   outputrehihihi has space allocated for dim+1 series of degree deg;
 *   outputrelohihi has space allocated for dim+1 series of degree deg;
 *   outputrehilohi has space allocated for dim+1 series of degree deg;
 *   outputrelolohi has space allocated for dim+1 series of degree deg;
 *   outputrehihilo has space allocated for dim+1 series of degree deg;
 *   outputrelohilo has space allocated for dim+1 series of degree deg;
 *   outputrehilolo has space allocated for dim+1 series of degree deg;
 *   outputrelololo has space allocated for dim+1 series of degree deg;
 *   outputimhihihi has space allocated for dim+1 series of degree deg;
 *   outputimlohihi has space allocated for dim+1 series of degree deg;
 *   outputimhilohi has space allocated for dim+1 series of degree deg;
 *   outputimlolohi has space allocated for dim+1 series of degree deg;
 *   outputimhihilo has space allocated for dim+1 series of degree deg;
 *   outputimlohilo has space allocated for dim+1 series of degree deg;
 *   outputimhilolo has space allocated for dim+1 series of degree deg;
 *   outputimlololo has space allocated for dim+1 series of degree deg;
 *   forwardrehihihi is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrelohihi is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrehilohi is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrelolohi is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardrehihilo is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrelohilo is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardrehilolo is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardrelololo is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimhihihi is work space for the highest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimlohihi is work space for the second highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimhilohi is work space for the third highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimlolohi is work space for the fourth highest doubles of nvr
 *                  forward  products, each product has deg+1 doubles;
 *   forwardimhihilo is work space for the fourth lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimlohilo is work space for the third lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   forwardimhilolo is work space for the second lowest doubles of nvr
 *                  forward products, each product deg+1 doubles;
 *   forwardimlololo is work space for the lowest doubles of nvr
 *                  forward products, each product has deg+1 doubles;
 *   backwardrehihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrehilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardrelololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhihihi is work space for the highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlohihi is work space for the second highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhilohi is work space for the third highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlolohi is work space for the fourth highest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhihilo is work space for the fourth lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlohilo is work space for the third lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimhilolo is work space for the second lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   backwardimlololo is work space for the lowest doubles of nvr-1
 *                  backward products, each product has deg+1 doubles;
 *   crossrehihihi  is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrelohihi  is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrehilohi  is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrelolohi  is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossrehihilo  is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrelohilo  is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrehilolo  is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossrelololo  is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimhihihi  is work space for the highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimlohihi  is work space for the second highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimhilohi  is work space for the third highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimlolohi  is work space for the fourth highest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles;
 *   crossimhihilo  is work space for the fourthlowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimlohilo  is work space for the third lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimhilolo  is work space for the second lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   crossimlololo  is work space for the lowest doubles of nvr-2
 *                  cross products, each product has deg+1 doubles.
 *   job            defines the convolution job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrehihihi are updated highest real parts of forward products;
 *   forwardrelohihi are updated second highest real parts of forward products;
 *   forwardrehilohi are updated third highest real parts of forward products;
 *   forwardrelolohi are updated fourth highest real parts of forward products;
 *   forwardrehihilo are updated fourth lowest real parts of forward products;
 *   forwardrelohilo are updated third lowest real parts of forward products;
 *   forwardrehilolo are updated second lowest real parts of forward products;
 *   forwardrelololo are updated lowest real parts forward products;
 *   forwardimhihihi are updated highest imag parts of forward products;
 *   forwardimlohihi are updated second highest imag parts of forward products;
 *   forwardimhilohi are updated third highest imag parts of forward products;
 *   forwardimlolohi are updated fourth highest imag parts of forward products;
 *   forwardimhihilo are updated fourth lowest imag parts of forward products;
 *   forwardimlohilo are updated third lowest imag parts of forward products;
 *   forwardimhilolo are updated second lowest imag parts of forward products;
 *   forwardimlololo are updated lowest imag parts forward products;
 *   backwardrehihihi are updated highest real parts of backward products;
 *   backwardrelohihi are updated 2nd highest real parts of backward products;
 *   backwardrehilohi are updated 3rd highest real parts of backward products;
 *   backwardrelolohi are updated 4th highest real parts of backward products;
 *   backwardrehihilo are updated 4th lowest real parts of backward products;
 *   backwardrelohilo are updated 3rd lowest real parts of backward products;
 *   backwardrehilolo are updated 2nd lowest real parts of backward products;
 *   backwardrelololo are updated lowest real parts backward products;
 *   backwardimhihihi are updated highest imag parts of backward products;
 *   backwardimlohihi are updated 2nd highest imag parts of backward products;
 *   backwardimhilohi are updated 3rd highest imag parts of backward products;
 *   backwardimlolohi are updated 4th highest imag parts of backward products;
 *   backwardimhihilo are updated 4th lowest imag parts of backward products;
 *   backwardimlohilo are updated 3rd lowest imag parts of backward products;
 *   backwardimhilolo are updated 2nd lowest imag parts of backward products;
 *   backwardimlololo are updated lowest parts backward products;
 *   crossrehihihi  are updated highest real parts of cross products;
 *   crossrelohihi  are updated second highest real parts of cross products;
 *   crossrehilohi  are updated third highest real parts of cross products;
 *   crossrelolohi  are updated fourth highest real parts of cross products;
 *   crossrehihilo  are updated fourth lowest real parts of cross products;
 *   crossrelohilo  are updated third lowest real parts of cross products;
 *   crossrehilolo  are updated second lowest real parts of cross products;
 *   crossrelololo  are updated lowest real parts cross products;
 *   crossimhihihi  are updated highest imag parts of cross products;
 *   crossimlohihi  are updated second highest imag parts of cross products;
 *   crossimhilohi  are updated third highest imag parts of cross products;
 *   crossimlolohi  are updated fourth highest imag parts of cross products;
 *   crossimhihilo  are updated fourth lowest imag parts of cross products;
 *   crossimlohilo  are updated third lowest imag parts of cross products;
 *   crossimhilolo  are updated second lowest imag parts of cross products;
 *   crossimlololo  are updated lowest imag parts cross products. */

void CPU_dbl8_add_job
 ( int deg,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job on real data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   csthihihi      highest parts of constant coefficient series;
 *   cstlohihi      second higest parts of constant coefficient series;
 *   csthilohi      third higest parts of constant coefficient series;
 *   cstlolohi      fourth higest parts of constant coefficient series;
 *   csthihilo      fourth lowest parts of constant coefficient series;
 *   cstlohilo      third lowest parts of constant coefficient series;
 *   csthilolo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cfflohihi      cffhilohi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhihilo[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cffhilolo[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cfflohihi[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflolohi[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cfflohilo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   forwardhihihi  are all highest parts of computed forward products;
 *   forwardlohihi  are all second highest parts of computed forward products;
 *   forwardhilohi  are all third highest parts of computed forward products;
 *   forwardlolohi  are all fourth highest parts of computed forward products;
 *   forwardhihilo  are all fourth lowest parts of computed forward products;
 *   forwardlohilo  are all third lowest parts of computed forward products;
 *   forwardhilolo  are all second lowest parts of computed forward products;
 *   forwardlololo  are all lowest parts of computed forward products;
 *   backwardhihihi are all highest parts of computed backward products;
 *   backwardlohihi are all second highest parts of computed backward products;
 *   backwardhilohi are all third highest parts of computed backward products;
 *   backwardlolohi are all fourth highest parts of computed backward products;
 *   backwardhihilo are all fourth lowest parts of computed backward products;
 *   backwardlohilo are all third lowest parts of computed backward products;
 *   backwardhilolo are all second lowest parts of computed backward products;
 *   backwardlololo are all lowest parts of computed backward products;
 *   crosshihihi    are all highest parts of computed cross products;
 *   crosslohihi    are all second highest parts of computed cross products;
 *   crosshilohi    are all third highest parts of computed cross products;
 *   crosslolohi    are all fourth highest parts of computed cross products;
 *   crosshihilo    are all fourth lowest parts of computed cross products;
 *   crosslohilo    are all third lowest parts of computed cross products;
 *   crosshilolo    are all second lowest parts of computed cross products;
 *   crosslololo    are all lowest parts of computed cross products;
 *   job            defines the addition job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihihi  are the updated highest parts of forward products;
 *   forwardlohihi  are the updated second highest parts of forward products;
 *   forwardhilohi  are the updated third highest parts of forward products;
 *   forwardlolohi  are the updated fourth highest parts of forward products;
 *   forwardhihilo  are the updated fourth lowest parts of forward products;
 *   forwardlohilo  are the updated third lowest parts of forward products;
 *   forwardhilolo  are the updated second lowest parts of forward products;
 *   forwardlololo  are the updated lowest parts forward products;
 *   backwardhihihi are the updated highest parts of backward products;
 *   backwardlohihi are the updated second highest parts of backward products;
 *   backwardhilohi are the updated third highest parts of backward products;
 *   backwardlolohi are the updated fourth highest parts of backward products;
 *   backwardhihilo are the updated fourth lowest parts of backward products;
 *   backwardlohilo are the updated third lowest parts of backward products;
 *   backwardhilolo are the updated second lowest parts of backward products;
 *   backwardlololo are the updated lowest parts backward products;
 *   crosshihihi    are the updated highest parts of cross products;
 *   crosslohihi    are the updated second highest parts of cross products;
 *   crosshilohi    are the updated third highest parts of cross products;
 *   crosslolohi    are the updated fourth highest parts of cross products;
 *   crosshihilo    are the updated fourth lowest parts of cross products;
 *   crosslohilo    are the updated third lowest parts of cross products;
 *   crosshilolo    are the updated second lowest parts of cross products;
 *   crosslololo    are the updated lowest parts cross products. */

void CPU_cmplx8_add_job
 ( int deg,
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
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job on complex data.
 *
 * ON ENTRY :
 *   deg            degree of the series;
 *   cstrehihihi    highest parts of real constant coefficients;
 *   cstrelohihi    second higest parts of real constant coefficients;
 *   cstrehilohi    third higest parts of real constant coefficients;
 *   cstrelolohi    fourth higest parts of real constant coefficients;
 *   cstrehihilo    fourth lowest parts of real constant coefficients;
 *   cstrelohilo    third lowest parts of real constant coefficients;
 *   cstrehilolo    second lowest parts of real constant coefficients;
 *   cstrelololo    lowest parts of real constant coefficients;
 *   cstimhihihi    highest parts of imag constant coefficients;
 *   cstimlohihi    second higest parts of imag constant coefficients;
 *   cstimhilohi    third higest parts of imag constant coefficients;
 *   cstimlolohi    fourth higest parts of imag constant coefficients;
 *   cstimhihilo    fourth lowest parts of imag constant coefficients;
 *   cstimlohilo    third lowest parts of imag constant coefficients;
 *   cstimhilolo    second lowest parts of imag constant coefficients;
 *   cstimlololo    lowest parts of imag constant coefficients;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrehilohi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilohi    cffrehihilo[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrehilolo[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrelohihi[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelolohi[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrelohilo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimhilohi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilohi    cffimhihilo[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimhilolo[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimlohihi[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlolohi[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimlohilo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   forwardrehihihi  are all highest real parts of forward products;
 *   forwardrelohihi  are all 2nd highest real parts of forward products;
 *   forwardrehilohi  are all 3rd highest real parts of forward products;
 *   forwardrelolohi  are all 4th highest real parts of forward products;
 *   forwardrehihilo  are all 4th lowest real parts of forward products;
 *   forwardrelohilo  are all 3rd lowest real parts of forward products;
 *   forwardrehilolo  are all 2nd lowest real parts of forward products;
 *   forwardrelololo  are all lowest real parts of forward products;
 *   forwardimhihihi  are all highest imag parts of forward products;
 *   forwardimlohihi  are all 2nd highest imag parts of forward products;
 *   forwardimhilohi  are all 3rd highest imag parts of forward products;
 *   forwardimlolohi  are all 4th highest imag parts of forward products;
 *   forwardimhihilo  are all 4th lowest imag parts of forward products;
 *   forwardimlohilo  are all 3rd lowest imag parts of forward products;
 *   forwardimhilolo  are all 2nd lowest imag parts of forward products;
 *   forwardimlololo  are all lowest imag parts of forward products;
 *   backwardrehihihi are all highest real parts of backward products;
 *   backwardrelohihi are all 2nd highest real parts of backward products;
 *   backwardrehilohi are all 3rd highest real parts of backward products;
 *   backwardrelolohi are all 4th highest real parts of backward products;
 *   backwardrehihilo are all 4th lowest real parts of backward products;
 *   backwardrelohilo are all 3rd lowest real parts of backward products;
 *   backwardrehilolo are all 2nd lowest real parts of backward products;
 *   backwardrelololo are all lowest real parts of backward products;
 *   backwardimhihihi are all highest imag parts of backward products;
 *   backwardimlohihi are all 2nd highest imag parts of backward products;
 *   backwardimhilohi are all 3rd highest imag parts of backward products;
 *   backwardimlolohi are all 4th highest imag parts of backward products;
 *   backwardimhihilo are all 4th lowest imag parts of backward products;
 *   backwardimlohilo are all 3rd lowest imag parts of backward products;
 *   backwardimhilolo are all 2nd lowest imag parts of backward products;
 *   backwardimlololo are all lowest imag parts of backward products;
 *   crossrehihihi  are all highest real parts of cross products;
 *   crossrelohihi  are all 2nd highest real parts of cross products;
 *   crossrehilohi  are all 3rd highest real parts of cross products;
 *   crossrelolohi  are all 4th highest real parts of cross products;
 *   crossrehihilo  are all 4th lowest real parts of cross products;
 *   crossrelohilo  are all 3rd lowest real parts of cross products;
 *   crossrehilolo  are all 2nd lowest real parts of cross products;
 *   crossrelololo  are all lowest real parts of cross products;
 *   crossimhihihi  are all highest imag parts of cross products;
 *   crossimlohihi  are all 2nd highest imag parts of cross products;
 *   crossimhilohi  are all 3rd highest imag parts of cross products;
 *   crossimlolohi  are all 4th highest imag parts of cross products;
 *   crossimhihilo  are all 4th lowest imag parts of cross products;
 *   crossimlohilo  are all 3rd lowest imag parts of cross products;
 *   crossimhilolo  are all 2nd lowest imag parts of cross products;
 *   crossimlololo  are all lowest imag parts of cross products;
 *   job            defines the addition job;
 *   verbose        if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrehihihi  are updated highest real parts of forward products;
 *   forwardrelohihi  are updated 2nd highest real parts of forward products;
 *   forwardrehilohi  are updated 3rd highest real parts of forward products;
 *   forwardrelolohi  are updated 4th highest real parts of forward products;
 *   forwardrehihilo  are updated 4th lowest real parts of forward products;
 *   forwardrelohilo  are updated 3rd lowest real parts of forward products;
 *   forwardrehilolo  are updated 2nd lowest real parts of forward products;
 *   forwardrelololo  are updated lowest real parts forward products;
 *   forwardimhihihi  are updated highest imag parts of forward products;
 *   forwardimlohihi  are updated 2nd highest imag parts of forward products;
 *   forwardimhilohi  are updated 3rd highest imag parts of forward products;
 *   forwardimlolohi  are updated 4th highest imag parts of forward products;
 *   forwardimhihilo  are updated 4th lowest imag parts of forward products;
 *   forwardimlohilo  are updated 3rd lowest imag parts of forward products;
 *   forwardimhilolo  are updated 2nd lowest imag parts of forward products;
 *   forwardimlololo  are updated lowest imag parts forward products;
 *   backwardrehihihi are updated highest real parts of backward products;
 *   backwardrelohihi are updated 2nd highest real parts of backward products;
 *   backwardrehilohi are updated 3rd highest real parts of backward products;
 *   backwardrelolohi are updated 4th highest real parts of backward products;
 *   backwardrehihilo are updated 4th lowest real parts of backward products;
 *   backwardrelohilo are updated 3rd lowest real parts of backward products;
 *   backwardrehilolo are updated 2nd lowest real parts of backward products;
 *   backwardrelololo are updated lowest real parts backward products;
 *   backwardimhihihi are updated highest imag parts of backward products;
 *   backwardimlohihi are updated 2nd highest imag parts of backward products;
 *   backwardimhilohi are updated 3rd highest imag parts of backward products;
 *   backwardimlolohi are updated 4th highest imag parts of backward products;
 *   backwardimhihilo are updated 4th lowest imag parts of backward products;
 *   backwardimlohilo are updated 3rd lowest imag parts of backward products;
 *   backwardimhilolo are updated 2nd imag lowest parts of backward products;
 *   backwardimlololo are updated lowest imag parts backward products;
 *   crossrehihihi  are updated highest real parts of cross products;
 *   crossrelohihi  are updated 2nd highest real parts of cross products;
 *   crossrehilohi  are updated 3rd highest real parts of cross products;
 *   crossrelolohi  are updated 4th highest real parts of cross products;
 *   crossrehihilo  are updated 4th lowest real parts of cross products;
 *   crossrelohilo  are updated 3rd lowest real parts of cross products;
 *   crossrehilolo  are updated 2nd lowest real parts of cross products;
 *   crossrelololo  are updated lowest real parts cross products;
 *   crossimhihihi  are updated highest imag parts of cross products;
 *   crossimlohihi  are updated 2nd highest imag parts of cross products;
 *   crossimhilohi  are updated 3rd highest imag parts of cross products;
 *   crossimlolohi  are updated 4th highest imag parts of cross products;
 *   crossimhihilo  are updated 4th lowest imag parts of cross products;
 *   crossimlohilo  are updated 3rd lowest imag parts of cross products;
 *   crossimhilolo  are updated 2nd lowest imag parts of cross products;
 *   crossimlololo  are updated lowest imag parts cross products. */

void CPU_dbl8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output
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
 *   csthihihi      highest parts of constant coefficient series;
 *   cstlohihi      second highest parts of constant coefficient series;
 *   csthilohi      third highest parts of constant coefficient series;
 *   cstlolohi      fourth highest parts of constant coefficient series;
 *   csthihilo      fourth lowest parts of constant coefficient series;
 *   cstlohilo      third lowest parts of constant coefficient series;
 *   csthilolo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest parts
 *                  of the coefficient series of monomial k;
 *   forwardhihihi  are all highest parts of forward products;
 *   forwardlohihi  are all second highest parts of forward products;
 *   forwardhilohi  are all third highest parts of forward products;
 *   forwardlolohi  are all fourth highest parts of forward products;
 *   forwardhihilo  are all fourth lowest parts of forward products;
 *   forwardlohilo  are all third lowest parts of forward products;
 *   forwardlohilo  are all second lowest parts of forward products;
 *   forwardlololo  are all lowest parts of forward products;
 *   backwardhihihi are all highest parts of backward products;
 *   backwardlohihi are all second highest parts of backward products;
 *   backwardhilohi are all third highest parts of backward products;
 *   backwardhilolo are all fourth highest parts of backward products;
 *   backwardhihilo are all fourth lowest parts of backward products;
 *   backwardlohilo are all third lowest parts of backward products;
 *   backwardhilolo are all second lowest parts of backward products;
 *   backwardlololo are all lowest parts of backward products;
 *   crosshihihi    are all highest parts of cross products;
 *   crossilohihi   are all second highest parts of cross products;
 *   crosshilohi    are all third highest parts of cross products;
 *   crosslolohi    are all fourth highest parts of cross products;
 *   crosshihilo    are all fourth lowest parts of cross products;
 *   crosslohilo    are all third lowest parts of cross products;
 *   crosshilolo    are all second lowest parts of cross products;
 *   crosslololo    are all lowest parts of cross products.
 *
 * ON RETURN :
 *   outputhihihi   highest parts of the values and all derivatives;
 *   outputlohihi   second highest parts of the values and all derivatives;
 *   outputhilohi   third highest parts of the values and all derivatives;
 *   outputlolohi   fourth highest parts of the values and all derivatives;
 *   outputhihilo   fourth lowest parts of the values and all derivatives;
 *   outputlohilo   third lowest parts of the values and all derivatives;
 *   outputhilolo   second lowest parts of the values and all derivatives;
 *   outputlololo   lowest parts of the values and all derivatives. */

void CPU_cmplx8_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo );
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
 *   cstrehihihi    highest real parts of constant coefficients;
 *   cstrelohihi    2nd highest real parts of constant coefficients;
 *   cstrehilohi    3rd highest real parts of constant coefficients;
 *   cstrelolohi    4th highest real parts of constant coefficients;
 *   cstrehihilo    4th lowest real parts of constant coefficients;
 *   cstrelohilo    3rd lowest real parts of constant coefficients;
 *   cstrehilolo    2nd lowest real parts of constant coefficients;
 *   cstrelololo    lowest real parts of constant coefficients;
 *   cstimhihihi    highest imag parts of constant coefficients;
 *   cstimlohihi    2nd highest imag parts of constant coefficients;
 *   cstimhilohi    3rd highest imag parts of constant coefficients;
 *   cstimlolohi    4th highest imag parts of constant coefficients;
 *   cstimhihilo    4th lowest imag parts of constant coefficients;
 *   cstimlohilo    3rd lowest imag parts of constant coefficients;
 *   cstimhilolo    2nd lowest imag parts of constant coefficients;
 *   cstimlololo    lowest imag parts of constant coefficients;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrelohihi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilohi    cffrehilohi[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrelolohi[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrehihilo[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelohilo[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrehilolo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimlohihi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilohi    cffimhilohi[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimlolohi[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimhihilo[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlohilo[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimhilolo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   forwardrehihihi are highest real parts of forward products;
 *   forwardrelohihi are 2nd highest real parts of forward products;
 *   forwardrehilohi are 3rd highest real parts of forward products;
 *   forwardrelolohi are 4th highest real parts of forward products;
 *   forwardrehihilo are 4th lowest real parts of forward products;
 *   forwardrelohilo are 3rd lowest real parts of forward products;
 *   forwardrelohilo are 2nd lowest real parts of forward products;
 *   forwardrelololo are lowest real parts of forward products;
 *   forwardimhihihi are highest imag parts of forward products;
 *   forwardimlohihi are 2nd highest imag parts of forward products;
 *   forwardimhilohi are 3rd highest imag parts of forward products;
 *   forwardimlolohi are 4th highest imag parts of forward products;
 *   forwardimhihilo are 4th lowest imag parts of forward products;
 *   forwardimlohilo are 3rd lowest imag parts of forward products;
 *   forwardimlohilo are 2nd lowest imag parts of forward products;
 *   forwardimlololo are lowest imag parts of forward products;
 *   backwardrehihihi are highest real parts of backward products;
 *   backwardrelohihi are 2nd highest real parts of backward products;
 *   backwardrehilohi are 3rd highest real parts of backward products;
 *   backwardrehilolo are 4th highest real parts of backward products;
 *   backwardrehihilo are 4th lowest real parts of backward products;
 *   backwardrelohilo are 3rd lowest real parts of backward products;
 *   backwardrehilolo are 2nd lowest real parts of backward products;
 *   backwardrelololo are lowest real parts of backward products;
 *   backwardimhihihi are highest imag parts of backward products;
 *   backwardimlohihi are 2nd highest imag parts of backward products;
 *   backwardimhilohi are 3rd highest imag parts of backward products;
 *   backwardimhilolo are 4th highest imag parts of backward products;
 *   backwardimhihilo are 4th lowest imag parts of backward products;
 *   backwardimlohilo are 3rd lowest imag parts of backward products;
 *   backwardimhilolo are 2nd lowest imag parts of backward products;
 *   backwardimlololo are lowest imag parts of backward products;
 *   crossrehihihi  are highest real parts of cross products;
 *   crossrelohihi  are 2nd highest real parts of cross products;
 *   crossrehilohi  are 3rd highest real parts of cross products;
 *   crossrelolohi  are 4th highest real parts of cross products;
 *   crossrehihilo  are 4th lowest real parts of cross products;
 *   crossrelohilo  are 3rd lowest real parts of cross products;
 *   crossrehilolo  are 2nd lowest real parts of cross products;
 *   crossrelololo  are lowest real parts of cross products;
 *   crossimhihihi  are highest imag parts of cross products;
 *   crossimlohihi  are 2nd highest imag parts of cross products;
 *   crossimhilohi  are 3rd highest imag parts of cross products;
 *   crossimlolohi  are 4th highest imag parts of cross products;
 *   crossimhihilo  are 4th lowest imag parts of cross products;
 *   crossimlohilo  are 3rd lowest imag parts of cross products;
 *   crossimhilolo  are 2nd lowest imag parts of cross products;
 *   crossimlololo  are lowest imag parts of cross products.
 *
 * ON RETURN :
 *   outputrehihihi are highest real parts of values and derivatives;
 *   outputrelohihi are 2nd highest real parts of values and derivatives;
 *   outputrehilohi are 3rd highest real parts of values and derivatives;
 *   outputrelolohi are 4th highest real parts of values and derivatives;
 *   outputrehihilo are 4th lowest real parts of values and derivatives;
 *   outputrelohilo are 3rd lowest real parts of values and derivatives;
 *   outputrehilolo are 2nd lowest real parts of values and derivatives;
 *   outputrelololo are lowest real parts of values and derivatives;
 *   outputimhihihi are highest imag parts of values and derivatives;
 *   outputimlohihi are 2nd highest imag parts of values and derivatives;
 *   outputimhilohi are 3rd highest imag parts of values and derivatives;
 *   outputimlolohi are 4th highest imag parts of values and derivatives;
 *   outputimhihilo are 4th lowest imag parts of values and derivatives;
 *   outputimlohilo are 3rd lowest imag parts of values and derivatives;
 *   outputimhilolo are 2nd lowest imag parts of values and derivatives;
 *   outputimlololo are lowest imag parts of values and derivatives. */

void CPU_dbl8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardhihihi, double ***forwardlohihi,
   double ***forwardhilohi, double ***forwardlolohi,
   double ***forwardhihilo, double ***forwardlohilo,
   double ***forwardhilolo, double ***forwardlololo,
   double ***backwardhihihi, double ***backwardlohihi,
   double ***backwardhilohi, double ***backwardlolohi, 
   double ***backwardhihilo, double ***backwardlohilo,
   double ***backwardhilolo, double ***backwardlololo, 
   double ***crosshihihi, double ***crosslohihi,
   double ***crosshilohi, double ***crosslolohi,
   double ***crosshihilo, double ***crosslohilo,
   double ***crosshilolo, double ***crosslololo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs, on real data.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   csthihihi      highest parts of constant coefficients;
 *   cstlohihi      second highest parts of constant coefficients;
 *   csthilohi      third highest parts of constant coefficients;
 *   cstlolohi      fourth highest parts of constant coefficients;
 *   csthihilo      fourth lowest parts of constant coefficients;
 *   cstlohilo      third lowest parts of constant coefficients;
 *   csthilolo      second lowest parts of constant coefficients;
 *   cstlololo      lowest parts of constant coefficients;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest
 *                  parts of the coefficient series of monomial k;
 *   cfflohihi      cfflohihi[k] has deg+1 doubles for the second highest
 *                  parts of the coefficient series of monomial k;
 *   cffhilohi      cffhilohi[k] has deg+1 doubles for the third highest
 *                  parts of the coefficient series of monomial k;
 *   cfflolohi      cfflolohi[k] has deg+1 doubles for the fourth highest
 *                  parts of the coefficient series of monomial k;
 *   cffhihilo      cffhihilo[k] has deg+1 doubles for the fourth lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflohilo      cfflohilo[k] has deg+1 doubles for the third lowest
 *                  parts of the coefficient series of monomial k;
 *   cffhilolo      cffhilolo[k] has deg+1 doubles for the second lowest
 *                  parts of the coefficient series of monomial k;
 *   cfflololo      cfflololo[k] has deg+1 doubles for the lowest
 *                  parts of the coefficient series of monomial k;
 *   forwardhihihi  are highest parts of forward products;
 *   forwardlohihi  are second highest parts of forward products;
 *   forwardhilohi  are third highest parts of forward products;
 *   forwardlolohi  are fourth highest parts of forward products;
 *   forwardhihilo  are fourth lowest parts of forward products;
 *   forwardlohilo  are third lowest parts of forward products;
 *   forwardlohilo  are second lowest parts of forward products;
 *   forwardlololo  are lowest parts of forward products;
 *   backwardhihihi are highest parts of backward products;
 *   backwardlohihi are second highest parts of backward products;
 *   backwardhilohi are third highest parts of backward products;
 *   backwardlolohi are fourth highest parts of backward products;
 *   backwardhihilo are fourth lowest parts of backward products;
 *   backwardlohilo are third lowest parts of backward products;
 *   backwardhilolo are second lowest parts of backward products;
 *   backwardlololo are lowest parts of backward products;
 *   jobs           defines the addition jobs;
 *   verbose        if true, then output is written.
 *
 * ON RETURN :
 *   outputhihihi   highest parts of values and derivatives;
 *   outputlohihi   second highest parts of values and derivatives;
 *   outputhilohi   third highest parts of values and derivatives;
 *   outputlolohi   fourth highest parts of values and derivatives;
 *   outputhihilo   fourth lowest parts of values and derivatives;
 *   outputlohilo   third lowest parts of values and derivatives;
 *   outputhilolo   second lowest parts of values and derivatives;
 *   outputlololo   lowest parts of values and derivatives. */

void CPU_cmplx8_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double ***forwardrehihihi, double ***forwardrelohihi,
   double ***forwardrehilohi, double ***forwardrelolohi,
   double ***forwardrehihilo, double ***forwardrelohilo,
   double ***forwardrehilolo, double ***forwardrelololo,
   double ***forwardimhihihi, double ***forwardimlohihi,
   double ***forwardimhilohi, double ***forwardimlolohi,
   double ***forwardimhihilo, double ***forwardimlohilo,
   double ***forwardimhilolo, double ***forwardimlololo,
   double ***backwardrehihihi, double ***backwardrelohihi,
   double ***backwardrehilohi, double ***backwardrelolohi, 
   double ***backwardrehihilo, double ***backwardrelohilo,
   double ***backwardrehilolo, double ***backwardrelololo, 
   double ***backwardimhihihi, double ***backwardimlohihi,
   double ***backwardimhilohi, double ***backwardimlolohi, 
   double ***backwardimhihilo, double ***backwardimlohilo,
   double ***backwardimhilolo, double ***backwardimlololo, 
   double ***crossrehihihi, double ***crossrelohihi,
   double ***crossrehilohi, double ***crossrelolohi,
   double ***crossrehihilo, double ***crossrelohilo,
   double ***crossrehilolo, double ***crossrelololo,
   double ***crossimhihihi, double ***crossimlohihi,
   double ***crossimhilohi, double ***crossimlolohi,
   double ***crossimhihilo, double ***crossimlohilo,
   double ***crossimhilolo, double ***crossimlololo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs, on real data.
 *
 * ON ENTRY :
 *   dim            total number of variables;
 *   nbr            number of monomials, excluding the constant term;
 *   deg            degree of the series;
 *   nvr            nvr[k] holds the number of variables in monomial k;
 *   idx            idx[k] has as many indices as the value of nvr[k],
 *                  idx[k][i] defines the place of the i-th variable,
 *                  with input values in input[idx[k][i]];
 *   cstrehihihi    highest real parts of constant coefficients;
 *   cstrelohihi    second highest real parts of constant coefficients;
 *   cstrehilohi    third highest real parts of constant coefficients;
 *   cstrelolohi    fourth highest real parts of constant coefficients;
 *   cstrehihilo    fourth lowest real parts of constant coefficients;
 *   cstrelohilo    third lowest real parts of constant coefficients;
 *   cstrehilolo    second lowest real parts of constant coefficients;
 *   cstrelololo    lowest real parts of constant coefficients;
 *   cstimhihihi    highest imag parts of constant coefficients;
 *   cstimlohihi    second highest imag parts of constant coefficients;
 *   cstimhilohi    third highest imag parts of constant coefficients;
 *   cstimlolohi    fourth highest imag parts of constant coefficients;
 *   cstimhihilo    fourth lowest imag parts of constant coefficients;
 *   cstimlohilo    third lowest imag parts of constant coefficients;
 *   cstimhilolo    second lowest imag parts of constant coefficients;
 *   cstimlololo    lowest imag parts of constant coefficients;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrelohihi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilohi    cffrehilohi[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrelolohi[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrehihilo[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelohilo[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrehilolo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimlohihi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilohi    cffimhilohi[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimlolohi[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimhihilo[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlohilo[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimhilolo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   forwardrehihihi are highest real parts of forward products;
 *   forwardrelohihi are second highest real parts of forward products;
 *   forwardrehilohi are third highest real parts of forward products;
 *   forwardrelolohi are fourth highest real parts of forward products;
 *   forwardrehihilo are fourth lowest real parts of forward products;
 *   forwardrelohilo are third lowest real parts of forward products;
 *   forwardrelohilo are second lowest real parts of forward products;
 *   forwardrelololo are lowest real parts of forward products;
 *   forwardimhihihi are highest imag parts of forward products;
 *   forwardimlohihi are second highest imag parts of forward products;
 *   forwardimhilohi are third highest imag parts of forward products;
 *   forwardimlolohi are fourth highest imag parts of forward products;
 *   forwardimhihilo are fourth lowest imag parts of forward products;
 *   forwardimlohilo are third lowest imag parts of forward products;
 *   forwardimlohilo are second lowest imag parts of forward products;
 *   forwardimlololo are lowest imag parts of forward products;
 *   backwardrehihihi are highest real parts of backward products;
 *   backwardrelohihi are second highest real parts of backward products;
 *   backwardrehilohi are third highest real parts of backward products;
 *   backwardrelolohi are fourth highest real parts of backward products;
 *   backwardrehihilo are fourth lowest real parts of backward products;
 *   backwardrelohilo are third lowest real parts of backward products;
 *   backwardrehilolo are second lowest real parts of backward products;
 *   backwardrelololo are lowest real parts of backward products;
 *   backwardimhihihi are highest imag parts of backward products;
 *   backwardimlohihi are second highest imag parts of backward products;
 *   backwardimhilohi are third highest imag parts of backward products;
 *   backwardimlolohi are fourth highest imag parts of backward products;
 *   backwardimhihilo are fourth lowest imag parts of backward products;
 *   backwardimlohilo are third lowest imag parts of backward products;
 *   backwardimhilolo are second lowest imag parts of backward products;
 *   backwardimlololo are lowest imag parts of backward products;
 *   jobs           defines the addition jobs;
 *   verbose        if true, then output is written.
 *
 * ON RETURN :
 *   outputrehihihi are highest real parts of values and derivatives;
 *   outputrelohihi are second highest real parts of values and derivatives;
 *   outputrehilohi are third highest real parts of values and derivatives;
 *   outputrelolohi are fourth highest real parts of values and derivatives;
 *   outputrehihilo are fourth lowest real parts of values and derivatives;
 *   outputrelohilo are third lowest real parts of values and derivatives;
 *   outputrehilolo are second lowest real parts of values and derivatives;
 *   outputrelololo are lowest real parts of values and derivatives;
 *   outputimhihihi are highest imag parts of values and derivatives;
 *   outputimlohihi are second highest imag parts of values and derivatives;
 *   outputimhilohi are third highest imag parts of values and derivatives;
 *   outputimlolohi are fourth highest imag parts of values and derivatives;
 *   outputimhihilo are fourth lowest imag parts of values and derivatives;
 *   outputimlohilo are third lowest imag parts of values and derivatives;
 *   outputimhilolo are second lowest imag parts of values and derivatives;
 *   outputimlololo are lowest imag parts of values and derivatives. */

void CPU_dbl8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   double *elapsedsec, int vrblvl );
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
 *   cstlohihi      second highest parts of constant coefficient series;
 *   csthilohi      third highest parts of constant coefficient series;
 *   cstlolohi      fourth highest parts of constant coefficient series;
 *   csthihilo      fourth lowest parts of constant coefficient series;
 *   cstlohilo      third lowest parts of constant coefficient series;
 *   csthilolo      second lowest parts of constant coefficient series;
 *   cstlololo      lowest parts of constant coefficient series;
 *   cffhihihi      cffhihihi[k] has deg+1 doubles for the highest parts
 *                  of the coefficient series of monomial k;
 *   cfflohihi      cffhilohi[k] has deg+1 doubles for the second highest
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
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg;
 *   cnvjobs        convolution jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputhihihi   has the highest parts of derivatives and the value,
 *                  outputhihihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihihi[dim] contains the value of the polynomial;
 *   outputlohihi   has the second highest parts of derivatives and the value,
 *                  outputlohihi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohihi[dim] contains the value of the polynomial;
 *   outputhilohi   has the third highest parts of derivatives and the value,
 *                  outputhilohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilohi[dim] contains the value of the polynomial;
 *   outputlolohi   has the fourth highest parts of derivatives and the value,
 *                  outputlolohi[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlolohi[dim] contains the value of the polynomial;
 *   outputhihilo   has the fourth lowest parts of derivatives and the value,
 *                  outputhihilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhihilo[dim] contains the value of the polynomial;
 *   outputlohilo   has the third lowest parts of derivatives and the value,
 *                  outputlohilo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlohilo[dim] contains the value of the polynomial;
 *   outputhilolo   has the second lowest parts of derivatives and the value,
 *                  outputhilolo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputhilolo[dim] contains the value of the polynomial;
 *   outputlololo   has the lowest parts of derivatives and the value,
 *                  outputlololo[k], for k from 0 to dim-1, contains the
 *                  derivative with respect to the variable k;
 *                  outputlololo[dim] contains the value of the polynomial;
 *   elapsedsec     is the elapsed time in seconds. */


void CPU_cmplx8_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
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
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
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
 *   cstrehihihi    highest real parts of constant coefficients;
 *   cstrelohihi    second highest real parts of constant coefficients;
 *   cstrehilohi    third highest real parts of constant coefficients;
 *   cstrelolohi    fourth highest real parts of constant coefficients;
 *   cstrehihilo    fourth lowest real parts of constant coefficients;
 *   cstrelohilo    third lowest real parts of constant coefficients;
 *   cstrehilolo    second lowest real parts of constant coefficients;
 *   cstrelololo    lowest real parts of constant coefficients;
 *   cstimhihihi    highest imag parts of constant coefficients;
 *   cstimlohihi    second highest imag parts of constant coefficients;
 *   cstimhilohi    third highest imag parts of constant coefficients;
 *   cstimlolohi    fourth highest imag parts of constant coefficients;
 *   cstimhihilo    fourth lowest imag parts of constant coefficients;
 *   cstimlohilo    third lowest imag parts of constant coefficients;
 *   cstimhilolo    second lowest imag parts of constant coefficients;
 *   cstimlololo    lowest imag parts of constant coefficients;
 *   cffrehihihi    cffrehihihi[k] has deg+1 doubles for the highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrehilohi[k] has deg+1 doubles for the second highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehihilo    cffrehihilo[k] has deg+1 doubles for the third highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrehilolo    cffrehilolo[k] has deg+1 doubles for the fourth highest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohihi    cffrelohihi[k] has deg+1 doubles for the fourth lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelolohi    cffrelolohi[k] has deg+1 doubles for the third lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelohilo    cffrelohilo[k] has deg+1 doubles for the second lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffrelololo    cffrelololo[k] has deg+1 doubles for the lowest
 *                  real parts of the coefficient series of monomial k;
 *   cffimhihihi    cffimhihihi[k] has deg+1 doubles for the highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimhilohi[k] has deg+1 doubles for the second highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhihilo    cffimhihilo[k] has deg+1 doubles for the third highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimhilolo    cffimhilolo[k] has deg+1 doubles for the fourth highest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohihi    cffimlohihi[k] has deg+1 doubles for the fourth lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlolohi    cffimlolohi[k] has deg+1 doubles for the third lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlohilo    cffimlohilo[k] has deg+1 doubles for the second lowest
 *                  imag parts of the coefficient series of monomial k;
 *   cffimlololo    cffimlololo[k] has deg+1 doubles for the lowest
 *                  imag parts of the coefficient series of monomial k;
 *   inputrehihihi  has the highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilohi  has the second highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehihilo  has the third highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrehilolo  has the fourth highest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohihi  has the fourth lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelolohi  has the third lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelohilo  has the second lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputrelololo  has the lowest real parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihihi  has the highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilohi  has the second highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhihilo  has the third highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimhilolo  has the fourth highest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohihi  has the fourth lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlolohi  has the third lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlohilo  has the second lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   inputimlololo  has the lowest imag parts of the power series
 *                  for all variables in the polynomial;
 *   outputrehihihi has space for dim+1 series of degree deg;
 *   outputrelohihi has space for dim+1 series of degree deg;
 *   outputrehilohi has space for dim+1 series of degree deg;
 *   outputrelolohi has space for dim+1 series of degree deg;
 *   outputrehihilo has space for dim+1 series of degree deg;
 *   outputrelohilo has space for dim+1 series of degree deg;
 *   outputrehilolo has space for dim+1 series of degree deg;
 *   outputrelololo has space for dim+1 series of degree deg;
 *   outputimhihihi has space for dim+1 series of degree deg;
 *   outputimlohihi has space for dim+1 series of degree deg;
 *   outputimhilohi has space for dim+1 series of degree deg;
 *   outputimlolohi has space for dim+1 series of degree deg;
 *   outputimhihilo has space for dim+1 series of degree deg;
 *   outputimlohilo has space for dim+1 series of degree deg;
 *   outputimhilolo has space for dim+1 series of degree deg;
 *   outputimlololo has space for dim+1 series of degree deg;
 *   cnvjobs        convolution jobs organized in layers;
 *   addjobs        addition jobs organized in layers;
 *   vrblvl         is the verbose level.
 *
 * ON RETURN :
 *   outputrehihihi has the highest real parts of derivatives and value,
 *                  outputrehihihi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrehihihi[dim] has the value of the polynomial;
 *   outputrelohihi has the 2nd highest real parts of derivatives and value,
 *                  outputrelohihi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrelohihi[dim] has the value of the polynomial;
 *   outputrehilohi has the 3rd highest real parts of derivatives and value,
 *                  outputrehilohi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrehilohi[dim] has the value of the polynomial;
 *   outputrelolohi has the 4th highest real parts of derivatives and value,
 *                  outputrelolohi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrelolohi[dim] has the value of the polynomial;
 *   outputrehihilo has the 4th lowest real parts of derivatives and value,
 *                  outputrehihilo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrehihilo[dim] has the value of the polynomial;
 *   outputrelohilo has the 3rd lowest real parts of derivatives and value,
 *                  outputrelohilo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrelohilo[dim] has the value of the polynomial;
 *   outputrehilolo has the 2nd lowest real parts of derivatives and value,
 *                  outputrehilolo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrehilolo[dim] has the value of the polynomial;
 *   outputrelololo has the lowest real parts of derivatives and value,
 *                  outputrelololo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputrelololo[dim] has the value of the polynomial;
 *   outputimhihihi has the highest imag parts of derivatives and value,
 *                  outputimhihihi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimhihihi[dim] has the value of the polynomial;
 *   outputimlohihi has the 2nd highest imag parts of derivatives and value,
 *                  outputimlohihi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimlohihi[dim] has the value of the polynomial;
 *   outputimhilohi has the 3rd highest imag parts of derivatives and value,
 *                  outputimhilohi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimhilohi[dim] has the value of the polynomial;
 *   outputimlolohi has the 4th highest imag parts of derivatives and value,
 *                  outputimlolohi[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimlolohi[dim] has the value of the polynomial;
 *   outputimhihilo has the 4th lowest imag parts of derivatives and value,
 *                  outputimhihilo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimhihilo[dim] has the value of the polynomial;
 *   outputimlohilo has the 3rd lowest imag parts of derivatives and value,
 *                  outputimlohilo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimlohilo[dim] has the value of the polynomial;
 *   outputimhilolo has the 2nd lowest imag parts of derivatives and value,
 *                  outputimhilolo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimhilolo[dim] has the value of the polynomial;
 *   outputimlololo has the lowest imag parts of derivatives and value,
 *                  outputimlololo[k], for k from 0 to dim-1, has the
 *                  derivative with respect to the variable k;
 *                  outputimlololo[dim] has the value of the polynomial;
 *   elapsedsec     is the elapsed time in seconds. */

#endif
