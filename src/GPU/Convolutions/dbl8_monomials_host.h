/* The file dbl8_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in octo double precision. */

#ifndef __dbl8_monomials_host_h__
#define __dbl8_monomials_host_h__

void CPU_dbl8_speel
 ( int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi, double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo, double *cffhilolo, double *cfflololo,
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
   double **crosshilolo, double **crosslololo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   multiplied with a coefficient series of the same degree,
 *   for real coefficients in octo double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr            number of variables in the product;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffhihihi      deg+1 doubles for the highest coefficient doubles;
 *   cfflohihi      deg+1 doubles for the second highest coefficient doubles;
 *   cffhilohi      deg+1 doubles for the third highest coefficient doubles;
 *   cfflolohi      deg+1 doubles for the fourth highest coefficient doubles;
 *   cffhihilo      deg+1 doubles for the fourth lowest coefficient doubles;
 *   cfflohilo      deg+1 doubles for the third lowest coefficient doubles;
 *   cffhilolo      deg+1 doubles for the second lowest coefficient doubles;
 *   cfflololo      deg+1 doubles for the lowest coefficient doubles;
 *   inputhihihi    stores the highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohihi    stores the second highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilohi    stores the third highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlolohi    stores the fourth highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhihilo    stores the fourth lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohilo    stores the third lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilolo    stores the second lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlololo    stores the lowest doubles of the input series
 *                  for all variables in the monomial;
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
 *   backwardhihihi is work space for the highest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohihi is work space for the second highest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilohi is work space for the third highest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardlolohi is work space for the fourth highest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardhihilo is work space for the fourth lowest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardlohilo is work space for the third lowest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardhilolo is work space for the second lowest doubles of nvr-2
 *                  backward products, each product has deg+1 doubles;
 *   backwardlololo is work space for the lowest doubles of nvr-2 
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
 *
 * ON RETURN :
 *   forwardhihihi  holds the highest doubles of the forward products,
 *   forwardlohihi  holds the second highest doubles of the forward products,
 *   forwardhilohi  holds the third highest doubles of the forward products,
 *   forwardlolohi  holds the fourth highest doubles of the forward products,
 *   forwardhihilo  holds the fourth lowest doubles of the forward products,
 *   forwardlohilo  holds the third lowest doubles of the forward products,
 *   forwardhilolo  holds the second lowest doubles of the forward products,
 *   forwardlololo  holds the lowest doubles of the forward products,
 *                  forward[nvr-1] contains the value of the product,
 *                  forward[nvr-2] contains the derivative with respect
 *                  to the last variable idx[nvr-1];
 *   backwardhihihi holds the highest doubles of the backward products,
 *   backwardlohihi holds the second highest doubles of the backward products,
 *   backwardhilohi holds the third highest doubles of the backward products,
 *   backwardlolohi holds the fourth highest doubles of the backward products,
 *   backwardhihilo holds the fourth lowest doubles of the backward products,
 *   backwardlohilo holds the third lowest doubles of the backward products,
 *   backwardhilolo holds the second lowest doubles of the backward products,
 *   backwardlololo holds the lowest doubles of the backward products,
 *                  backward[nvr-3] contains the derivative with respect
 *                  to the first variable idx[0];
 *   crosshihihi    holds the highest doubles of the cross products,
 *   crosslohihi    holds the second highest doubles of the cross products,
 *   crosshilohi    holds the third highest doubles of the cross products,
 *   crosslolohi    holds the fourth highest doubles of the cross products,
 *   crosshihilo    holds the fourth lowest doubles of the cross products,
 *   crosslohilo    holds the third lowest doubles of the cross products,
 *   crosshilolo    holds the second lowest doubles of the cross products,
 *   crosslololo    holds the lowest doubles of the cross products,
 *                  cross[k] contains the derivatve with respect to
 *                  variable idx[k+1]. */

void CPU_cmplx8_speel
 ( int nvr, int deg, int *idx,
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
   double **crossimhilolo, double **crossimlololo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in octo double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr              number of variables in the product;
 *   deg              truncation degree of the series;
 *   idx              as many indices as the value of nvr,
 *                    idx[k] defines the place of the k-th variable,
 *                    with input values in input[idx[k]];
 *   cffrehihihi      highest doubles of the real parts of the coefficients
 *                    of the series of the product;
 *   cffrelohihi      second highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehilohi      third highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelolohi      fourth highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehihilo      fourth lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelohilo      third lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehilolo      second lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelololo      lowest doubles of the real parts of the coefficients
 *                    of the series of the product;
 *   cffimhihihi      highest doubles of the imaginary parts of the 
 *                    coefficients of the series of the product;
 *   cffimlohihi      second highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhilohi      third highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlolohi      fourth highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhihilo      fourth lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlohilo      third lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhilolo      second lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlololo      lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
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
 *   inputimhilolo    holds the third lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimhilolo    holds the second lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   forwardrehihihi  is work space for the highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrelohihi  is work space for the second highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrehilohi  is work space for the third highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrelolohi  is work space for the fourth highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrehihilo  is work space for the fourth lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrelohilo  is work space for the third lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrehilolo  is work space for the second lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardrelololo  is work space for the lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimhihihi  is work space for the highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimlohihi  is work space for the second highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimhilohi  is work space for the third highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimlolohi  is work space for the fourth highest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimhihilo  is work space for the fourth lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimlohilo  is work space for the third lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimhilolo  is work space for the second lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   forwardimlololo  is work space for the lowest doubles of nvr
 *                    forward products, each product has deg+1 doubles;
 *   backwardrehihihi is work space for the highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrelohihi is work space for the second highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrehilohi is work space for the third highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrelolohi is work space for the fourth highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrehihilo is work space for the fourth lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrelohilo is work space for the third lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrehilolo is work space for the second lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardrelololo is work space for the lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimhihihi is work space for the highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimlohihi is work space for the second highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimhilohi is work space for the third highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimlolohi is work space for the fourth highest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimhihilo is work space for the fourth lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimlohilo is work space for the third lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimhilolo is work space for the second lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   backwardimlololo is work space for the lowest doubles of nvr-2
 *                    backward products, each product has deg+1 doubles;
 *   crossrehihihi    is work space for the highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrelohihi    is work space for the second highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrehilohi    is work space for the third highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrelolohi    is work space for the fourth highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrehihilo    is work space for the fourth lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrelohilo    is work space for the third lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrehilolo    is work space for the second lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossrelololo    is work space for the lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossimhihihi    is work space for the highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossimlohihi    is work space for the second highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossimhilohi    is work space for the third highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossimlolohi    is work space for the fourth highest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles;
 *   crossimhihilo    is work space for the fourth lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles.
 *   crossimlohilo    is work space for the third lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles.
 *   crossimhilolo    is work space for the second lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles.
 *   crossimlololo    is work space for the lowest doubles of nvr-2
 *                    cross products, each product has deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrehihihi  holds the highest doubles of the real parts
 *                    of the forward products,
 *   forwardrelohihi  holds the second highest doubles of the real parts
 *                    of the forward products,
 *   forwardrehilohi  holds the third highest doubles of the real parts
 *                    of the forward products,
 *   forwardrelolohi  holds the fourth highest doubles of the real parts
 *                    of the forward products,
 *   forwardrehihilo  holds the fourth lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrelohilo  holds the third lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrehilolo  holds the second lowest doubles of the real parts
 *                    of the forward products,
 *   forwardrelololo  holds the lowest doubles of the real parts
 *                    of the forward products,
 *   forwardimhihihi  holds the highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlohihi  holds the second highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhilohi  holds the third highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlolohi  holds the fourth highest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhihilo  holds the fourth lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlohilo  holds the third lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimhilolo  holds the second lowest doubles of the imaginary parts
 *                    of the forward products,
 *   forwardimlololo  holds the lowest doubles of the imaginary parts
 *                    of the forward products,
 *                    forward[nvr-1] contains the value of the product,
 *                    forward[nvr-2] contains the derivative with respect
 *                    to the last variable idx[nvr-1];
 *   backwardrehihihi holds the highest doubles of the real parts
 *                    of the backward products,
 *   backwardrelohihi holds the second highest doubles of the real parts
 *                    of the backward products,
 *   backwardrehilohi holds the third highest doubles of the real parts
 *                    of the backward products,
 *   backwardrelolohi holds the fourth highest doubles of the real parts
 *                    of the backward products,
 *   backwardrehihilo holds the fourth lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrelohilo holds the third lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrehilolo holds the second lowest doubles of the real parts
 *                    of the backward products,
 *   backwardrelololo holds the lowest doubles of the real parts
 *                    of the backward products,
 *   backwardimhihihi holds the highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimlohihi holds the second highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimhilohi holds the third highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimlolohi holds the fourth highest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimhihilo holds the fourth lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimlohilo holds the third lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimhilolo holds the second lowest doubles of the imaginary parts
 *                    of the backward products,
 *   backwardimlololo holds the lowest doubles of the imaginary parts
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
 *                    cross[k] contains the derivative with respect to
 *                    the variable idx[k+1]. */

void CPU_dbl8_evaldiff
 ( int dim, int nvr, int deg, int *idx,
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
 *   dim            total number of variables;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffhihihi      deg+1 doubles for the highest coefficient doubles;
 *   cfflohihi      deg+1 doubles for the second highest coefficient doubles;
 *   cffhilohi      deg+1 doubles for the third highest coefficient doubles;
 *   cfflolohi      deg+1 doubles for the fourth highest coefficient doubles;
 *   cffhihilo      deg+1 doubles for the fourth lowest coefficient doubles;
 *   cfflohilo      deg+1 doubles for the third lowest coefficient doubles;
 *   cffhilolo      deg+1 doubles for the second lowest coefficient doubles;
 *   cfflololo      deg+1 doubles for the lowest coefficient doubles;
 *   inputhihihi    stores the highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohihi    stores the second highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilohi    stores the third highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlolohi    stores the fourth highest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhihilo    stores the fourth lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlohilo    stores the third lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputhilolo    stores the second lowest doubles of the input series
 *                  for all variables in the monomial;
 *   inputlololo    stores the lowest doubles of the input series
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
 *   outputhihihi   highest doubles of the derivatives and the value,
 *   outputlohihi   second highest doubles of the derivatives and the value,
 *   outputhilohi   third highest doubles of the derivatives and the value,
 *   outputlolohi   fourth highest doubles of the derivatives and the value,
 *   outputhihilo   fourth lowest doubles of the derivatives and the value,
 *   outputlohilo   third lowest doubles of the derivatives and the value,
 *   outputhilolo   second lowest doubles of the derivatives and the value,
 *   outputlololo   lowest doubles of the derivatives and the value,
 *                  output[idx[k]], for k from 0 to nvr, contains the
 *                  deriviative with respect to the variable idx[k];
 *                  output[dim] contains the value of the product. */

void CPU_cmplx8_evaldiff
 ( int dim, int nvr, int deg, int *idx,
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
 *   dim              total number of variables;
 *   deg              truncation degree of the series;
 *   idx              as many indices as the value of nvr,
 *                    idx[k] defines the place of the k-th variable,
 *                    with input values in input[idx[k]];
 *   cffrehihihi      highest doubles of the real parts of the coefficients
 *                    of the series of the product;
 *   cffrelohihi      second highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehilohi      third highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelolohi      fourth highest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehihilo      fourth lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelohilo      third lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrehilolo      second lowest doubles of the real parts of the
 *                    coefficients of the series of the product;
 *   cffrelololo      lowest doubles of the real parts of the coefficients
 *                    of the series of the product;
 *   cffimhihihi      highest doubles of the imaginary parts of the 
 *                    coefficients of the series of the product;
 *   cffimlohihi      second highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhilohi      third highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlolohi      fourth highest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhihilo      fourth lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlohilo      third lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimhilolo      second lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
 *   cffimlololo      lowest doubles of the imaginary parts of the
 *                    coefficients of the series of the product;
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
 *   inputimhilolo    holds the third lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimhilolo    holds the second lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   inputimlololo    holds the lowest doubles of the imaginary parts
 *                    of the input series for all variables;
 *   outputrehihihi   has space allocated for dim+1 series of degree deg;
 *   outputrelohihi   has space allocated for dim+1 series of degree deg;
 *   outputrehilohi   has space allocated for dim+1 series of degree deg;
 *   outputrelolohi   has space allocated for dim+1 series of degree deg;
 *   outputrehihilo   has space allocated for dim+1 series of degree deg;
 *   outputrelohilo   has space allocated for dim+1 series of degree deg;
 *   outputrehilolo   has space allocated for dim+1 series of degree deg;
 *   outputrelololo   has space allocated for dim+1 series of degree deg;
 *   outputimhihihi   has space allocated for dim+1 series of degree deg;
 *   outputimlohihi   has space allocated for dim+1 series of degree deg;
 *   outputimhilohi   has space allocated for dim+1 series of degree deg;
 *   outputimlolohi   has space allocated for dim+1 series of degree deg;
 *   outputimhihilo   has space allocated for dim+1 series of degree deg;
 *   outputimlohilo   has space allocated for dim+1 series of degree deg;
 *   outputimhilolo   has space allocated for dim+1 series of degree deg;
 *   outputimlololo   has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputrehihihi   stores the highest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrelohihi   stores the second highest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrehilohi   stores the third highest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrelolohi   stores the fourth highest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrehihilo   stores the fourth lowest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrelohilo   stores the third lowest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrehilolo   stores the second lowest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputrelololo   stores the lowest doubles of the real parts
 *                    of the derivatives and the value,
 *   outputimhihihi   stores the highest doubles of the imaginary parts
 *                    of the derivatives and the value,
 *   outputimlohihi   stores the second highest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimhilohi   stores the third highest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimlolohi   stores the fourth highest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimhihilo   stores the fourth lowest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimlohilo   stores the third lowest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimhilolo   stores the second lowest doubles of the imaginary
 *                    parts of the derivatives and the value,
 *   outputimlololo   stores the lowest doubles of the imaginary parts
 *                    of the derivatives and the value,
 *                    output[idx[k]], for k from 0 to nvr, stores the
 *                    deriviative with respect to the variable idx[k];
 *                    output[dim] stores the value of the product. */

#endif
