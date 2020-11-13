/* The file dbl4_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in quad double precision. */

#ifndef __dbl4_monomials_host_h__
#define __dbl4_monomials_host_h__

void CPU_dbl4_speel
 ( int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   multiplied with a coefficient series of the same degree,
 *   for real coefficients in quad double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr          number of variables in the product;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffhihi      deg+1 doubles for the highest doubles in the coefficients;
 *   cfflohi      deg+1 doubles for the second highest coefficient doubles;
 *   cffhilo      deg+1 doubles for the second lowest coefficient doubles;
 *   cfflolo      deg+1 doubles for the lowest doubles in the coefficients;
 *   inputhihi    contains the highest doubles of the input series
 *                for all variables in the monomial;
 *   inputlohi    contains the second highest doubles of the input series
 *                for all variables in the monomial;
 *   inputhilo    contains the second lowest doubles of the input series
 *                for all variables in the monomial;
 *   inputlolo    contains the lowest doubles of the input series
 *                for all variables in the monomial;
 *   forwardhihi  is work space for the highest doubles of nvr forward 
 *                products, forwardhihi[k] can store deg+1 doubles;
 *   forwardlohi  is work space for the second highest doubles of nvr
 *                forward  products, forwardlohi[k] can store deg+1 doubles;
 *   forwardhilo  iss work space for the second lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   forwardlolo  is work space for the lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   backwardhihi is work space for the highest doubles of nvr-2 backward
 *                products, backwardhihi[k] can store deg+1 doubles;
 *   backwardlohi is work space for the second highest doubles of nvr-2
 *                backward products, backwardlohi[k] can store deg+1 doubles;
 *   backwardhilo is work space for the second lowest doubles of nvr-2
 *                backward products, backwardhilo[k] can store deg+1 doubles;
 *   backwardlolo is work space for the lowest doubles of nvr-2 backward
 *                products, backwardlolo[k] can store deg+1 doubles;
 *   crosshihi    is work space for the highest doubles of nvr-2 cross
 *                products, crosshi[k] can store deg+1 doubles;
 *   crosslohi    is work space for the second highest doubles of nvr-2
 *                cross products, crosshi[k] can store deg+1 doubles;
 *   crosshilo    is work space for the second lowest doubles of nvr-2
 *                cross products, crosshilo[k] can store for deg+1 doubles.
 *   crosslolo    is work space for the lowest doubles of nvr-2 cross
 *                products, crosslolo[k] can store for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardhihi  holds the highest doubles of the forward products,
 *   forwardlohi  holds the second highest doubles of the forward products,
 *   forwardhilo  holds the second lowest doubles of the forward products,
 *   forwardlolo  holds the lowest doubles of the forward products,
 *                forward[nvr-1] contains the value of the product,
 *                forward[nvr-2] contains the derivative with respect
 *                to the last variable idx[nvr-1];
 *   backwardhihi holds the highest doubles of the backward products,
 *   backwardlohi holds the second highest doubles of the backward products,
 *   backwardhilo holds the second lowest doubles of the backward products,
 *   backwardlolo holds the lowest doubles of the backward products,
 *                backward[nvr-3] contains the derivative with respect
 *                to the first variable idx[0];
 *   crosshihi    holds the highest doubles of the cross products,
 *   crosslohi    holds the second highest doubles of the cross products,
 *   crosshilo    holds the second lowest doubles of the cross products,
 *   crosslolo    holds the lowest doubles of the cross products,
 *                cross[k] contains the derivatve with respect to
 *                variable idx[k+1]. */

void CPU_cmplx4_speel
 ( int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
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
   double **crossimhilo, double **crossimlolo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in quad double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr            number of variables in the product;
 *   deg            truncation degree of the series;
 *   idx            as many indices as the value of nvr,
 *                  idx[k] defines the place of the k-th variable,
 *                  with input values in input[idx[k]];
 *   cffrehihi      highest doubles of the real parts of the coefficients
 *                  of the series of the product;
 *   cffrelohi      second highest doubles of the real parts of the
 *                  coefficients of the series of the product;
 *   cffrehilo      second lowest doubles of the real parts of the
 *                  coefficients of the series of the product;
 *   cffrelolo      lowest doubles of the real parts of the coefficients
 *                  of the series of the product;
 *   cffimhihi      highest doubles of the imaginary parts of the coefficients
 *                  of the series of the product;
 *   cffimlohi      second highest doubles of the imaginary parts of the
 *                  coefficients of the series of the product;
 *   cffimhilo      second lowest doubles of the imaginary parts of the
 *                  coefficients of the series of the product;
 *   cffimlolo      lowest doubles of the imaginary parts of the coefficients
 *                  of the series of the product;
 *   inputrehihi    holds the highest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrelohi    holds the second highest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrehilo    holds the second lowest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrelolo    holds the lowest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputimhihi    holds the highest doubles of the imaginary parts of the
 *                  coefficients of the series for all variables;
 *   inputimlohi    holds the highest doubles of the imaginary parts of the
 *                  coefficients of the series for all variables;
 *   inputimhilo    holds the second lowest doubles of the imaginary parts of
 *                  the coefficients of the series for all variables;
 *   inputimlolo    holds the lowest doubles of the imaginary parts of the
 *                  coefficients of the series for all variables;
 *   forwardrehihi  is work space for the highest doubles of nvr forward
 *                  products, forwardrehihi[k] holds deg+1 doubles;
 *   forwardrelohi  is work space for the second highest doubles of nvr
 *                  forward products, forwardrehihi[k] holds deg+1 doubles;
 *   forwardrehilo  is work space for the second lowest doubles of nvr forward
 *                  products, forwardrehilo[k] holds deg+1 doubles;
 *   forwardrelolo  is work space for the lowest doubles of nvr forward
 *                  products, forwardrelolo[k] holds deg+1 doubles;
 *   forwardimhihi  is work space for the highest doubles of nvr forward
 *                  products, forwardimhihi[k] holds deg+1 doubles;
 *   forwardimlohi  is work space for the second highest doubles of nvr
 *                  forward products, forwardimlohi[k] holds deg+1 doubles;
 *   forwardimhilo  is work space for the second lowest doubles of nvr forward
 *                  products, forwardimlo[k] holds deg+1 doubles;
 *   forwardimlolo  is work space for the lowest doubles of nvr forward
 *                  products, forwardimlolo[k] holds deg+1 doubles;
 *   backwardrehihi is work space for the highest doubles of nvr-2 backward
 *                  products, backwardrehihi[k] holds deg+1 doubles;
 *   backwardrelohi is work space for the second highest doubles of nvr-2
 *                  backward products, backwardrelohi[k] holds deg+1 doubles;
 *   backwardrehilo is work space for the second lowest doubles of nvr-2
 *                  backward products, backwardrehilo[k] holds deg+1 doubles;
 *   backwardrelolo is work space for the lowest doubles of nvr-2 backward
 *                  products, backwardrelolo[k] holds deg+1 doubles;
 *   backwardimhihi is work space for the highest doubles of nvr-2 backward
 *                  products, backwardimhihi[k] holds deg+1 doubles;
 *   backwardimlohi is work space for the second highest doubles of nvr-2
 *                  backward products, backwardimlohi[k] holds deg+1 doubles;
 *   backwardimhilo is work space for the second lowest doubles of nvr-2
 *                  backward products, backwardimhilo[k] holds deg+1 doubles;
 *   backwardimlolo is work space for the lowest doubles of nvr-2 backward
 *                  products, backwardimlolo[k] holds deg+1 doubles;
 *   crossrehihi    is work space for the highest doubles of nvr-2 cross
 *                  products, crossrehihi[k] holds deg+1 doubles;
 *   crossrelohi    is work space for the second highest doubles of nvr-2
 *                  cross products, crossrelohi[k] holds deg+1 doubles;
 *   crossrehilo    is work space for the second lowest doubles of nvr-2 cross
 *                  products, crossrehilo[k] holds deg+1 doubles;
 *   crossrelolo    is work space for the lowest doubles of nvr-2 cross
 *                  products, crossrelolo[k] holds deg+1 doubles;
 *   crossimhihi    is work space for the highest doubles of nvr-2 cross
 *                  products, crossimhihi[k] holds deg+1 doubles;
 *   crossimlohi    is work space for the second highest doubles of nvr-2
 *                  cross products, crossimlohi[k] holds deg+1 doubles;
 *   crossimhilo    is work space for the second lowest doubles of nvr-2 cross
 *                  products, crossimlolo[k] holds deg+1 doubles.
 *   crossimlolo    is work space for the lowest doubles of nvr-2 cross
 *                  products, crossimlolo[k] holds deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrehihi  holds the highest doubles of the real parts
 *                  of the forward products,
 *   forwardrelohi  holds the second highest doubles of the real parts
 *                  of the forward products,
 *   forwardrehilo  holds the second lowest doubles of the real parts
 *                  of the forward products,
 *   forwardrelolo  holds the lowest doubles of the real parts
 *                  of the forward products,
 *   forwardimhihi  holds the highest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimlohi  holds the second highest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimhilo  holds the second lowest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimlolo  holds the lowest doubles of the imaginary parts
 *                  of the forward products,
 *                  forward[nvr-1] contains the value of the product,
 *                  forward[nvr-2] contains the derivative with respect
 *                  to the last variable idx[nvr-1];
 *   backwardrehihi holds the highest doubles of the real parts
 *                  of the backward products,
 *   backwardrelohi holds the second highest doubles of the real parts
 *                  of the backward products,
 *   backwardrehilo holds the second lowest doubles of the real parts
 *                  of the backward products,
 *   backwardrelolo holds the lowest doubles of the real parts
 *                  of the backward products,
 *   backwardimhihi holds the highest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardimlohi holds the second highest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardimhilo holds the second lowest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardimlolo holds the lowest doubles of the imaginary parts
 *                  of the backward products,
 *                  backward[nvr-3] contains the derivative with respect
 *                  to the first variable idx[0];
 *   crossrehihi    stores the highest doubles of the real parts
 *                  of the cross products,
 *   crossrelohi    stores the second highest doubles of the real parts
 *                  of the cross products,
 *   crossrehilo    stores the second lowest doubles of the real parts
 *                  of the cross products,
 *   crossrelolo    stores the lowest doubles of the real parts
 *                  of the cross products,
 *   crossimhihi    stores the highest doubles of the imaginary parts
 *                  of the cross products,
 *   crossimlohi    stores the second highest doubles of the imaginary parts
 *                  of the cross products,
 *   crossimhilo    stores the second lowest doubles of the imaginary parts
 *                  of the cross products,
 *   crossimlolo    stores the lowest doubles of the imaginary parts
 *                  of the cross products,
 *                  cross[k] contains the derivative with respect to
 *                  the variable idx[k+1]. */

void CPU_dbl4_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffhihi      deg+1 doubles for the highest doubles in the coefficients;
 *   cfflohi      deg+1 doubles for the second highest coefficient doubles;
 *   cffhilo      deg+1 doubles for the second lowest coefficient doubles;
 *   cfflolo      deg+1 doubles for the lowest doubles in the coefficients;
 *   inputhihi    contains the highest doubles of the input series
 *                for all variables in the monomial;
 *   inputlohi    contains the second highest doubles of the input series
 *                for all variables in the monomial;
 *   inputhilo    contains the second lowest doubles of the input series
 *                for all variables in the monomial;
 *   inputlolo    contains the lowest doubles of the input series
 *                for all variables in the monomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputhihi   holds the highest doubles of the derivatives and the value,
 *   outputlohi   second highest doubles of the derivatives and the value,
 *   outputhilo   second lowest doubles of the derivatives and the value,
 *   outputlolo   holds the lowest doubles of the derivatives and the value,
 *                output[idx[k]], for k from 0 to nvr, contains the
 *                deriviative with respect to the variable idx[k];
 *                output[dim] contains the value of the product. */

void CPU_cmplx4_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo );
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
 *   cffrehihi      highest doubles of the real parts of the coefficients
 *                  of the series of the product;
 *   cffrelohi      second highest doubles of the real parts of the
 *                  coefficients of the series of the product;
 *   cffrehilo      second lowest doubles of the real parts of the
 *                  coefficients of the series of the product;
 *   cffrelolo      lowest doubles of the real parts of the coefficients
 *                  of the series of the product;
 *   cffimhihi      highest doubles of the imaginary parts of the coefficient
 *                  of the series of the product;
 *   cffimlohi      second highest doubles of the imaginary parts of the
 *                  coefficient of the series of the product;
 *   cffimhilo      second lowest doubles of the imaginary parts of the
 *                  coefficient of the series of the product;
 *   cffimlolo      lowest doubles of the imaginary parts of the coefficient
 *                  of the series of the product;
 *   inputrehihi    holds the highest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrelohi    holds the second highest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrehilo    holds the second lowest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputrelolo    holds the lowest doubles of the real parts of the
 *                  coefficients of the series for all variables;
 *   inputimhihi    holds the highest doubles of the imaginary parts of the
 *                  coefficients of the series for all variables;
 *   inputimlohi    holds the second highest doubles of the imaginary parts
 *                  of the coefficients of the series for all variables;
 *   inputimhilo    holds the second lowest doubles of the imaginary parts of
 *                  the coefficients of the series for all variables;
 *   inputimlolo    holds the lowest doubles of the imaginary parts of the
 *                  coefficients of the series for all variables;
 *   outputrehihi   has space allocated for dim+1 series of degree deg;
 *   outputrelohi   has space allocated for dim+1 series of degree deg;
 *   outputrehilo   has space allocated for dim+1 series of degree deg;
 *   outputrelolo   has space allocated for dim+1 series of degree deg;
 *   outputimhihi   has space allocated for dim+1 series of degree deg;
 *   outputimlohi   has space allocated for dim+1 series of degree deg;
 *   outputimhilo   has space allocated for dim+1 series of degree deg;
 *   outputimlolo   has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputrehihi   stores the highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelohi   stores the second highest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrehilo   stores the second lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputrelolo   stores the lowest doubles of the real parts
 *                  of the derivatives and the value,
 *   outputimhihi   stores the highest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *   outputimlohi   stores the second highest doubles of the imaginary
 *                  parts of the derivatives and the value,
 *   outputimhilo   stores the second lowest doubles of the imaginary
 *                  parts of the derivatives and the value,
 *   outputimlolo   stores the lowest doubles of the imaginary parts
 *                  of the derivatives and the value,
 *                  output[idx[k]], for k from 0 to nvr, stores the
 *                  deriviative with respect to the variable idx[k];
 *                  output[dim] stores the value of the product. */

#endif
