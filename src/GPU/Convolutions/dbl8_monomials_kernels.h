// The file dbl8_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in octo double precision.

#ifndef __dbl8_monomials_kernels_h__
#define __dbl8_monomials_kernels_h__

void GPU_dbl8_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo,
   double *inputhihihi, double *inputlohihi,
   double *inputhilohi, double *inputlolohi,
   double *inputhihilo, double *inputlohilo,
   double *inputhilolo, double *inputlololo,
   double *forwardhihihi, double *forwardlohihi,
   double *forwardhilohi, double *forwardlolohi,
   double *forwardhihilo, double *forwardlohilo,
   double *forwardhilolo, double *forwardlololo,
   double *backwardhihihi, double *backwardlohihi,
   double *backwardhilohi, double *backwardlolohi,
   double *backwardhihilo, double *backwardlohilo,
   double *backwardhilolo, double *backwardlololo,
   double *crosshihihi, double *crosslohihi,
   double *crosshilohi, double *crosslolohi,
   double *crosshihilo, double *crosslohilo,
   double *crosshilolo, double *crosslololo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in octo double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS             number of threads in one block, must be deg+1;
 *   nvr            number of variables in the product;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffhihihi      deg+1 highest doubles of the coefficient series;
 *   cfflohihi      deg+1 second highest doubles of the coefficient series;
 *   cffhilohi      deg+1 third highest doubles of the coefficient series;
 *   cfflolohi      deg+1 fourth highest doubles of the coefficient series;
 *   cffhihilo      deg+1 fourth lowest doubles of the coefficient series;
 *   cfflohilo      deg+1 third lowest doubles of the coefficient series;
 *   cffhilolo      deg+1 second lowest doubles of the coefficient series;
 *   cfflololo      deg+1 lowest doubles of the coefficient series;
 *   inputhihihi    holds the highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohihi    holds the second highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilohi    holds the third highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlolohi    holds the fourth highest of the input series
 *                  for all variables in the monomial;
 *   inputhihilo    holds the fourth lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohilo    holds the third lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilolo    holds the second lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlololo    holds the lowest doubles of the input series
 *                  for all variables in the monomial;
 *   forwardhihihi  is work space for the highest doubles of all nvr
 *                  forward products, forwardhihihi[k] holds deg+1 doubles;
 *   forwardlohihi  is work space for the second highest doubles of all nvr 
 *                  forward products, forwardlohihi[k] holds deg+1 doubles;
 *   forwardhilohi  is work space for the third highest doubles of all nvr 
 *                  forward products, forwardhilohi[k] holds deg+1 doubles;
 *   forwardlolohi  is work space for the fourth highest doubles of all nvr 
 *                  forward products, forwardlolohi[k] holds deg+1 doubles;
 *   forwardhihilo  is work space for the fourth lowest doubles of all nvr
 *                  forward products, forwardhihilo[k] holds deg+1 doubles;
 *   forwardlohilo  is work space for the third lowest doubles of all nvr
 *                  forward products, forwardhilolo[k] holds deg+1 doubles;
 *   forwardhilolo  is work space for the second lowest doubles of all nvr
 *                  forward products, forwardhilolo[k] holds deg+1 doubles;
 *   forwardlololo  is work space for the lowest doubles of all nvr
 *                  forward products, forwardlololo[k] holds deg+1 doubles;
 *   backwardhihihi is work space for the highest doubles of all nvr-2
 *                  backward products, backwardhihihi[k] holds deg+1 doubles;
 *   backwardlohihi is work space for the second highest doubles of all nvr-2
 *                  backward products, backwardlohihi[k] holds deg+1 doubles;
 *   backwardhilohi is work space for the third highest doubles of all nvr-2
 *                  backward products, backwardhilohi[k] holds deg+1 doubles;
 *   backwardlolohi is work space for the fourth highest doubles of all nvr-2
 *                  backward products, backwardlolohi[k] holds deg+1 doubles;
 *   backwardhihilo is work space for the fourth lowest doubles of all nvr-2
 *                  backward products, backwardhihilo[k] holds deg+1 doubles;
 *   backwardlohilo is work space for the third lowest doubles of all nvr-2
 *                  backward products, backwardlohilo[k] holds deg+1 doubles;
 *   backwardhilolo is work space for the second lowest doubles of all nvr-2
 *                  backward products, backwardhilolo[k] holds deg+1 doubles;
 *   backwardlololo is work space for the lowest doubles of all nvr-2
 *                  backward products, backwardlololo[k] holds deg+1 doubles;
 *   crosshihihi    is work space for the highest doubles of all nvr-2
 *                  cross products, crosshihihi[k] can hold deg+1 doubles;
 *   crosslohihi    is work space for the second highest doubles of all nvr-2
 *                  cross products, crosslohihi[k] can hold deg+1 doubles;
 *   crosshilohi    is work space for the third highest doubles of all nvr-2
 *                  cross products, crosshilohi[k] can hold deg+1 doubles;
 *   crosslolohi    is work space for the fourth highest doubles of all nvr-2
 *                  cross products, crosslolohi[k] can hold deg+1 doubles;
 *   crosshihilo    is work space for the fourth lowest doubles of all nvr-2
 *                  cross products, crosshihilo[k] can hold deg+1 doubles.
 *   crosslohilo    is work space for the third lowest doubles of all nvr-2
 *                  cross products, crosslohilo[k] can hold deg+1 doubles.
 *   crosshilolo    is work space for the second lowest doubles of all nvr-2
 *                  cross products, crosshilolo[k] can hold deg+1 doubles.
 *   crosslololo    is work space for the lowest doubles of all nvr-2 cross
 *                  products, crosslololo[k] can hold deg+1 doubles.
 *
 * ON RETURN :
 *   forwardhihihi  stores the highest doubles of the forward products,
 *   forwardlohihi  stores the second highest doubles of the forward products,
 *   forwardhilohi  stores the third highest doubles of the forward products,
 *   forwardlolohi  stores the fourth highest doubles of the forward products,
 *   forwardhihilo  stores the fourth lowest doubles of the forward products,
 *   forwardlohilo  stores the third lowest doubles of the forward products,
 *   forwardhilolo  stores the second lowest doubles of the forward products,
 *   forwardlololo  stores the lowest doubles of the forward products,
 *                  forward[nvr-1] contains the value of the product,
 *                  forward[nvr-2] contains the derivative with respect
 *                  to the last variable idx[nvr-1] if nvr > 2;
 *   backwardhihihi stores the highest doubles of the backward products,
 *   backwardlohihi stores the second highest doubles of the backward products,
 *   backwardhilohi stores the third highest doubles of the backward products,
 *   backwardlolohi stores the fourth highest doubles of the backward products,
 *   backwardhihilo stores the fourth lowest doubles of the backward products,
 *   backwardlohilo stores the third lowest doubles of the backward products,
 *   backwardhilolo stores the second lowest doubles of the backward products,
 *   backwardlololo stores the lowest doubles of the backward products,
 *                  backward[nvr-3] contains the derivative with respect
 *                  to the first variable idx[0] if nvr > 2;
 *   crosshihihi    stores the highest doubles of the cross products,
 *   crosslohihi    stores the second highest doubles of the cross products,
 *   crosshilohi    stores the third highest doubles of the cross products,
 *   crosslolohi    stores the fourth highest doubles of the cross products,
 *   crosshihilo    stores the fourth lowest doubles of the cross products,
 *   crosslohilo    stores the third lowest doubles of the cross products,
 *   crosshilolo    stores the second lowest doubles of the cross products,
 *   crosslololo    stores the lowest doubles of the cross products,
 *                  cross[k] contains the derivatve with respect to
 *                  variable idx[k+1]. */

void GPU_cmplx8_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double *inputrehihihi, double *inputrelohihi,
   double *inputrehilohi, double *inputrelolohi,
   double *inputrehihilo, double *inputrelohilo,
   double *inputrehilolo, double *inputrelololo,
   double *inputimhihihi, double *inputimlohihi,
   double *inputimhilohi, double *inputimlolohi,
   double *inputimhihilo, double *inputimlohilo,
   double *inputimhilolo, double *inputimlololo,
   double *forwardrehihihi, double *forwardrelohihi,
   double *forwardrehilohi, double *forwardrelolohi,
   double *forwardrehihilo, double *forwardrelohilo,
   double *forwardrehilolo, double *forwardrelololo,
   double *forwardimhihihi, double *forwardimlohihi,
   double *forwardimhilohi, double *forwardimlolohi,
   double *forwardimhihilo, double *forwardimlohilo,
   double *forwardimhilolo, double *forwardimlololo,
   double *backwardrehihihi, double *backwardrelohihi,
   double *backwardrehilohi, double *backwardrelolohi,
   double *backwardrehihilo, double *backwardrelohilo,
   double *backwardrehilolo, double *backwardrelololo,
   double *backwardimhihihi, double *backwardimlohihi,
   double *backwardimhilohi, double *backwardimlolohi,
   double *backwardimhihilo, double *backwardimlohilo,
   double *backwardimhilolo, double *backwardimlololo,
   double *crossrehihihi, double *crossrelohihi,
   double *crossrehilohi, double *crossrelolohi,
   double *crossrehihilo, double *crossrelohilo,
   double *crossrehilolo, double *crossrelololo,
   double *crossimhihihi, double *crossimlohihi,
   double *crossimhilohi, double *crossimlolohi,
   double *crossimhihilo, double *crossimlohilo,
   double *crossimhilolo, double *crossimlololo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in octo double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS               number of threads in one block, must be deg+1;
 *   nvr              number of variables in the product;
 *   deg              truncation degree of the series;
 *   idx              as many indices as the value of nvr,
 *                    idx[k] defines the place of the k-th variable,
 *                    with input values in input[idx[k]];
 *   cffrehihihi      highest doubles of real parts of coefficients;
 *   cffrelohihi      second highest doubles of real parts of coefficients;
 *   cffrehilohi      third highest doubles of real parts of coefficients;
 *   cffrelolohi      fourth highest doubles of real parts of coefficients;
 *   cffrehihilo      fourth lowest doubles of real parts of coefficients;
 *   cffrelohilo      third lowest doubles of real parts of coefficients;
 *   cffrehilolo      second lowest doubles of real parts of coefficients;
 *   cffrelololo      lowest doubles of the real parts of the coefficients;
 *   cffimhihihi      highest doubles of imaginary parts of the coefficients;
 *   cffimlohihi      second highest doubles of imaginary coefficients;
 *   cffimhilohi      third highest doubles of imaginary coefficients;
 *   cffimlolohi      fourth highest doubles of imaginary coefficients;
 *   cffimhihilo      fourth lowest doubles of imaginary coefficients;
 *   cffimlohilo      thir lowest doubles of imaginary coefficients;
 *   cffimhilolo      second lowest doubles of imaginary coefficients;
 *   cffimlololo      lowest doubles of imaginary parts of coefficients;
 *   inputrehihihi    holds the highest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrelohihi    holds the second highest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrehilohi    holds the third highest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrelolohi    holds the fourth highest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrehihilo    holds the fourth lowest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrelohilo    holds the third lowest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrehilolo    holds the second lowest doubles of the real parts
 *                    of the input series for all variables;
 *   inputrelololo    holds the lowest doubles of the real parts
 *                    of the input series for all variables;
 *   inputimhihihi    holds the highest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimlohihi    holds the second highest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimhilohi    holds the third highest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimlolohi    holds the fourth highest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimhihilo    holds the fourth lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimlohilo    holds the third lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimhilolo    holds the second lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimlololo    holds the lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   forwardrehihihi  is work space for the highest doubles of nvr forward
 *                    products, for all real parts of the coefficients,
 *                    forwardrehihi[k] has space for deg+1 doubles;
 *   forwardrelohihi  is work space for the second highest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrelohihi[k] has space for deg+1 doubles;
 *   forwardrehilohi  is work space for the third highest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrehilohi[k] has space for deg+1 doubles;
 *   forwardrelolohi  is work space for the fourth highest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrelolohi[k] has space for deg+1 doubles;
 *   forwardrehihilo  is work space for the fourth lowest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrehihilo[k] has space for deg+1 doubles;
 *   forwardrelohilo  is work space for the third lowest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrelohilo[k] has space for deg+1 doubles;
 *   forwardrehilolo  is work space for the second lowest doubles of nvr
 *                    forward products, for all real parts of the coefficients,
 *                    forwardrehilolo[k] has space for deg+1 doubles;
 *   forwardrelololo  is work space for the lowest doubles of nvr forward
 *                    products, for all real parts of the coefficients,
 *                    forwardrelololo[k] has space for deg+1 doubles;
 *   forwardimhihihi  is work space for the highest doubles of nvr forward
 *                    products, for all imaginary parts of the coefficients,
 *                    forwardimhihihi[k] has space for deg+1 doubles;
 *   forwardimlohihi  is work space for the second highest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimlohihi[k] has space for deg+1 doubles;
 *   forwardimhilohi  is work space for the third highest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimhilohi[k] has space for deg+1 doubles;
 *   forwardimlolohi  is work space for the fourth highest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimlolohi[k] has space for deg+1 doubles;
 *   forwardimhihilo  is work space for the fourth lowest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimhihilo[k] has space for deg+1 doubles;
 *   forwardimlohilo  is work space for the third lowest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimlohilo[k] has space for deg+1 doubles;
 *   forwardimhilolo  is work space for the second lowest doubles of nvr
 *                    forward products, for imaginary parts of coefficients,
 *                    forwardimhilolo[k] has space for deg+1 doubles;
 *   forwardimlololo  is work space for the lowest doubles of nvr forward
 *                    products, for all imaginary parts of the coefficients,
 *                    forwardimlololo[k] has space for deg+1 doubles;
 *   backwardrehihihi is work space for all highest doubles of nvr-2 backward
 *                    products, for all real parts of the coefficients,
 *                    backwardrehihihi[k] has space for deg+1 doubles;
 *   backwardrelohihi is work space for all second highest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrelohihi[k] has space for deg+1 doubles;
 *   backwardrehilohi is work space for all third highest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrehilohi[k] has space for deg+1 doubles;
 *   backwardrelolohi is work space for all fourth highest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrelolohi[k] has space for deg+1 doubles;
 *   backwardrehihilo is work space for all fourth lowest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrehihilo[k] has space for deg+1 doubles;
 *   backwardrelohilo is work space for all third lowest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrelohilo[k] has space for deg+1 doubles;
 *   backwardrehilolo is work space for all second lowest doubles of nvr-2
 *                    backward products, for all real parts of the coefficients,
 *                    backwardrehilolo[k] has space for deg+1 doubles;
 *   backwardrelololo is work space for all lowest doubles of nvr-2 backward
 *                    products, for all real parts of the coefficients,
 *                    backwardrelololo[k] has space for deg+1 doubles;
 *   backwardimhihihi is work space for the highest doubles of nvr-2 backward
 *                    products, for all imaginary parts of the coefficients,
 *                    backwardimhihihi[k] has space for deg+1 doubles;
 *   backwardimlohihi is work space for the second highest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimlohihi[k] has space for deg+1 doubles;
 *   backwardimhilohi is work space for the third highest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimhilohi[k] has space for deg+1 doubles;
 *   backwardimlolohi is work space for the fourth highest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimlolohi[k] has space for deg+1 doubles;
 *   backwardimhihilo is work space for the fourth lowest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimhihilo[k] has space for deg+1 doubles;
 *   backwardimlohilo is work space for the third lowest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimhilo[k] has space for deg+1 doubles;
 *   backwardimhilolo is work space for the second lowest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimhilolo[k] has space for deg+1 doubles;
 *   backwardimlololo is work space for the lowest doubles of nvr-2
 *                    backward products, for imaginary parts of coefficients,
 *                    backwardimlololo[k] has space for deg+1 doubles;
 *   crossrehihihi    is work space for the highest doubles of nvr-2 cross
 *                    products, for the real parts of the coefficients,
 *                    crossrehihihi[k] has space for deg+1 doubles;
 *   crossrelohihi    is work space for the second highest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrelohihi[k] has space for deg+1 doubles;
 *   crossrehilohi    is work space for the third highest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrehilohi[k] has space for deg+1 doubles;
 *   crossrelolohi    is work space for the fourth highest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrelolohi[k] has space for deg+1 doubles;
 *   crossrehihilo    is work space for the fourth lowest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrehihilo[k] has space for deg+1 doubles;
 *   crossrelohilo    is work space for the second lowest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrelohilo[k] has space for deg+1 doubles;
 *   crossrehilolo    is work space for the second lowest doubles of nvr-2
 *                    cross products, for the real parts of the coefficients,
 *                    crossrehilolo[k] has space for deg+1 doubles;
 *   crossrelololo    is work space for the lowest doubles of nvr-2 cross
 *                    products, for the real parts of the coefficients,
 *                    crossrelololo[k] has space for deg+1 doubles;
 *   crossimhihihi    is work space for the highest doubles of nvr-2 cross
 *                    products, for the imaginary parts of the coefficients,
 *                    crossimhihihi[k] has space for deg+1 doubles.
 *   crossimlohihi    is work space for the second highest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimlohihi[k] has space for deg+1 doubles.
 *   crossimhilohi    is work space for the third highest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimhilohi[k] has space for deg+1 doubles.
 *   crossimlolohi    is work space for the fourth highest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimlolohi[k] has space for deg+1 doubles.
 *   crossimhihilo    is work space for the fourth lowest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimhihilo[k] has space for deg+1 doubles;
 *   crossimlohilo    is work space for the third lowest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimlohilo[k] has space for deg+1 doubles;
 *   crossimhilolo    is work space for the second lowest doubles of nvr-2
 *                    cross products, for the imaginary parts of coefficients,
 *                    crossimhilolo[k] has space for deg+1 doubles;
 *   crossimlololo    is work space for the lowest doubles of nvr-2 cross
 *                    products, for the imaginary parts of the coefficients,
 *                    crossimlololo[k] has space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrehihihi  stores the highest doubles of the real parts
 *                    of the forward products,
 *   forwardrelohihi  stores the second highest doubles of the real parts
 *                    of the forward products,
 *   forwardrehilohi  stores the third highest doubles of the real parts
 *                    of the forward products,
 *   forwardrelolohi  stores the fourth highest doubles of the real parts
 *                    of the forward products,
 *   forwardrehihilo  stores the fourth lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrelohilo  stores the third lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrehilolo  stores the second lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrelololo  stores the lowest doubles of the real parts
 *                    of the forward products,
 *   forwardimhihihi  stores the highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlohihi  stores the second highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhilohi  stores the third highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlolohi  stores the fourth highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhihilo  stores the fourth lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlohilo  stores the third lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhilolo  stores the second lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlololo  stores the lowest doubles of the imaginary parts
 *                    of the forward products,
 *                    forward[nvr-1] contains the value of the product,
 *                    forward[nvr-2] contains the derivative with respect
 *                    to the last variable idx[nvr-1];
 *   backwardrehihihi stores the highest doubles of the real parts
 *                    of the backward products,
 *   backwardrelohihi stores the second highest doubles of the real parts
 *                    of the backward products,
 *   backwardrehilohi stores the third highest doubles of the real parts
 *                    of the backward products,
 *   backwardrelolohi stores the fourth highest doubles of the real parts
 *                    of the backward products,
 *   backwardrehihilo stores the fourth lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrelohilo stores the third lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrehilolo stores the second lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrelololo stores the lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrehihihi stores the highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrelohihi stores the second highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrehilohi stores the third highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrelolohi stores the fourth highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrehihilo stores the third lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrelohilo stores the second lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrehilolo stores the second lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardrelololo stores the lowest doubles of the imaginary parts
 *                    of the backward products,
 *                    backward[nvr-3] contains the derivative with respect
 *                    to the first variable idx[0];
 *   crossrehihihi    stores the highest doubles of the real parts
 *                    of the cross products,
 *   crossrelohihi    stores the second highest doubles of the real parts
 *                    of the cross products,
 *   crossrehilohi    stores the third highest doubles of the real parts
 *                    of the cross products,
 *   crossrelolohi    stores the fourth highest doubles of the real parts
 *                    of the cross products,
 *   crossrehihilo    stores the fourth lowest doubles of the real parts
 *                    of the cross products,
 *   crossrelohilo    stores the third lowest doubles of the real parts
 *                    of the cross products,
 *   crossrehilolo    stores the second lowest doubles of the real parts
 *                    of the cross products,
 *   crossrelololo    stores the lowest doubles of the real parts
 *                    of the cross products,
 *   crossimhihihi    stores the highest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimlohihi    stores the second highest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimhilohi    stores the third highest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimlolohi    stores the fourth highest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimhihilo    stores the fourth lowest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimlohilo    stores the third lowest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimhilolo    stores the second lowest doubles of the imaginary parts
 *                    of the cross products,
 *   crossimlololo    stores the lowest doubles of the imaginary parts
 *                    of the cross products,
 *                    cross[k] contains the derivatve with respect to
 *                    variable idx[k+1]. */

void GPU_dbl8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi, double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo, double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   BS             number of threads in one block, must be deg+1;
 *   dim            total number of variables;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffhihihi      deg+1 highest doubles of the coefficient series;
 *   cfflohihi      deg+1 second highest doubles of the coefficient series;
 *   cffhilohi      deg+1 third highest doubles of the coefficient series;
 *   cfflolohi      deg+1 fourth highest doubles of the coefficient series;
 *   cffhihilo      deg+1 fourth lowest doubles of the coefficient series;
 *   cfflohilo      deg+1 third lowest doubles of the coefficient series;
 *   cffhilolo      deg+1 second lowest doubles of the coefficient series;
 *   cfflololo      deg+1 lowest doubles of the coefficient series;
 *   inputhihihi    holds the highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohihi    holds the second highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilohi    holds the third highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlolohi    holds the fourth highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhihilo    holds the fourth lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohilo    holds the third lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilolo    holds the second lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlololo    holds the lowest doubles of the input series
 *                  for all variables in the monomial;
 *   outputhihihi   has space allocated for dim+1 series of degree deg;
 *   outputlohihi   has space allocated for dim+1 series of degree deg;
 *   outputhilohi   has space allocated for dim+1 series of degree deg;
 *   outputlolohi   has space allocated for dim+1 series of degree deg;
 *   outputhihilo   has space allocated for dim+1 series of degree deg;
 *   outputlohilo   has space allocated for dim+1 series of degree deg;
 *   outputhilolo   has space allocated for dim+1 series of degree deg;
 *   outputlololo   has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputhihihi   stores the highest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputlohihi   stores the second highest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputhilohi   stores the third highest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputlolohi   stores the fourth highest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputhihilo   stores the fourth lowest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputlohilo   stores the third lowest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputhilolo   stores the second lowest doubles of the derivatives and 
 *                  the value of the monomial,
 *   outputlololo   stores the lowest doubles of the derivatives and 
 *                  the value of the monomial,
 *                  output[idx[k]], for k from 0 to nvr, contains the
 *                  deriviative with respect to the variable idx[k];
 *                  output[dim] contains the value of the product. */

void GPU_cmplx8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   BS             number of threads in one block, must be deg+1;
 *   dim            total number of variables;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffrehihihi    highest doubles of real parts of coefficients;
 *   cffrelohihi    second highest doubles of real parts of coefficients;
 *   cffrehilohi    third highest doubles of real parts of coefficients;
 *   cffrelolohi    fourth highest doubles of real parts of coefficients;
 *   cffrehihilo    fourth lowest doubles of real parts of coefficients;
 *   cffrelohilo    third lowest doubles of real parts of coefficients;
 *   cffrehilolo    second lowest doubles of real parts of coefficients;
 *   cffrelololo    lowest doubles of the real parts of the coefficients;
 *   cffimhihihi    highest doubles of imaginary parts of the coefficients;
 *   cffimlohihi    second highest doubles of imaginary coefficients;
 *   cffimhilohi    third highest doubles of imaginary coefficients;
 *   cffimlolohi    fourth highest doubles of imaginary coefficients;
 *   cffimhihilo    fourth lowest doubles of imaginary coefficients;
 *   cffimlohilo    third lowest doubles of imaginary coefficients;
 *   cffimhilolo    second lowest doubles of imaginary coefficients;
 *   cffimlololo    lowest doubles of imaginary parts of coefficients;
 *   inputrehihihi  holds the highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelohihi  holds the second highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrehilohi  holds the third highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelolohi  holds the fourth highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrehihilo  holds the fourth lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelohilo  holds the third lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrehilolo  holds the second lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelololo  holds the lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputimhihihi  holds the highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlohihi  holds the second highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimhilohi  holds the third highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlolohi  holds the fourth highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimhihilo  holds the fourth lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlohilo  holds the third lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimhilolo  holds the second lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlololo  holds the lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   outputrehihihi has space allocated for dim+1 series of degree deg, 
 *                  for the highest doubles of the real output;
 *   outputrelohihi has space allocated for dim+1 series of degree deg,
 *                  for the third second highest doubles of the real output;
 *   outputrehilohi has space allocated for dim+1 series of degree deg,
 *                  for the third highest doubles of the real output;
 *   outputrelolohi has space allocated for dim+1 series of degree deg,
 *                  for the fourth highest doubles of the real output;
 *   outputrehihilo has space allocated for dim+1 series of degree deg,
 *                  for the fourth lowest doubles of the real output;
 *   outputrelohilo has space allocated for dim+1 series of degree deg,
 *                  for the third lowest doubles of the real output;
 *   outputrehilolo has space allocated for dim+1 series of degree deg,
 *                  for the second lowest doubles of the real output;
 *   outputrelololo has space allocated for dim+1 series of degree deg,
 *                  for the lowest doubles of the real utput;
 *   outputimhihihi has space allocated for dim+1 series of degree deg,
 *                  for the highest doubles of the imaginary output;
 *   outputimlohihi has space allocated for dim+1 series of degree deg,
 *                  for the second highest doubles of the imaginary output;
 *   outputimhilohi has space allocated for dim+1 series of degree deg,
 *                  for the third highest doubles of the imaginary output;
 *   outputimlolohi has space allocated for dim+1 series of degree deg,
 *                  for the fourth highest doubles of the imaginary output;
 *   outputimhihilo has space allocated for dim+1 series of degree deg, 
 *                  for the fourth lowest doubles of the imaginary output;
 *   outputimlohilo has space allocated for dim+1 series of degree deg, 
 *                  for the third lowest doubles of the imaginary output;
 *   outputimhilolo has space allocated for dim+1 series of degree deg, 
 *                  for the second lowest doubles of the imaginary output;
 *   outputimlololo has space allocated for dim+1 series of degree deg,
 *                  for the lowest doubles of the imaginary output.
 *
 * ON RETURN :
 *   outputrehihihi stores the highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelohihi stores the second highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrehilohi stores the third highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelolohi stores the fourth highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrehihilo stores the fourth lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelohilo stores the third lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrehilolo stores the second lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelololo stores the lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputimhihihi stores the highest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimlohihi stores the second highest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimhilohi stores the third highest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimlolohi stores the fourth highest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimhihilo stores the fourth lowest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimlohilo stores the third lowest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimhilolo stores the second lowest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimlololo stores the lowest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *                  output[idx[k]], for k from 0 to nvr, contains the
 *                  derivative with respect to the variable idx[k];
 *                  output[dim] contains the value of the product. */

#endif
