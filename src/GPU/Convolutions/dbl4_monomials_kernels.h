// The file dbl4_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in quad double precision.

#ifndef __dbl4_monomials_kernels_h__
#define __dbl4_monomials_kernels_h__

void GPU_dbl4_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double *inputhihi, double *inputlohi, double *inputhilo, double *inputlolo,
   double *forwardhihi, double *forwardlohi,
   double *forwardhilo, double *forwardlolo,
   double *backwardhihi, double *backwardlohi,
   double *backwardhilo, double *backwardlolo,
   double *crosshihi, double *crosslohi,
   double *crosshilo, double *crosslolo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in quad double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS           number of threads in one block, must be deg+1;
 *   nvr          number of variables in the product;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffhihi      deg+1 highest doubles of the coefficient series;
 *   cfflohi      deg+1 second highest doubles of the coefficient series;
 *   cffhilo      deg+1 second lowest doubles of the coefficient series;
 *   cfflolo      deg+1 lowest doubles of the coefficient series;
 *   inputhihi    holds the highest doubles of the input series
 *                for all variables in the monomial;
 *   inputlohi    holds the second highest doubles of the input series
 *                for all variables in the monomial;
 *   inputhilo    holds the second lowest doubles of the input series
 *                for all variables in the monomial;
 *   inputlolo    holds the lowest doubles of the input series
 *                for all variables in the monomial;
 *   forwardhihi  has work space for the highest doubles of all nvr forward
 *                products, forwardhihi[k] can hold deg+1 doubles;
 *   forwardlohi  has work space for the second highest doubles of all nvr 
 *                forward products, forwardhi[k] can hold deg+1 doubles;
 *   forwardhilo  has work space for the second lowest doubles of all nvr
 *                forward products, forwardhilo[k] can hold deg+1 doubles;
 *   forwardlolo  has work space for the lowest doubles of all nvr forward
 *                products, forwardlolo[k] can hold deg+1 doubles;
 *   backwardhihi has work space for the highest doubles of all nvr-2 backward
 *                products, backwardhihi[k] can hold deg+1 doubles;
 *   backwardlohi has work space for the second highest doubles of all nvr-2
 *                backward products, backwardlohi[k] can hold deg+1 doubles;
 *   backwardhilo has work space for the second lowest doubles of all nvr-2
 *                backward products, backwardhilo[k] can hold deg+1 doubles;
 *   backwardlolo has work space for the lowest doubles of all nvr-2 backward
 *                products, backwardlolo[k] can hold deg+1 doubles;
 *   crosshihi    has work space for the highest doubles of all nvr-2 cross
 *                products, crosshihi[k] can hold deg+1 doubles;
 *   crosslohi    has work space for the second highest doubles of all nvr-2
 *                cross products, crosslohi[k] can hold deg+1 doubles;
 *   crosshilo    has work space for the second lowest doubles of all nvr-2
 *                cross products, crosshilo[k] can hold deg+1 doubles.
 *   crosslolo    has work space for the lowest doubles of all nvr-2 cross
 *                products, crosslolo[k] can hold deg+1 doubles.
 *
 * ON RETURN :
 *   forwardhihi  stores the highest doubles of the forward products,
 *   forwardlohi  stores the second highest doubles of the forward products,
 *   forwardhilo  stores the second lowest doubles of the forward products,
 *   forwardlolo  stores the lowest doubles of the forward products,
 *                forward[nvr-1] contains the value of the product,
 *                forward[nvr-2] contains the derivative with respect
 *                to the last variable idx[nvr-1] if nvr > 2;
 *   backwardhihi stores the highest doubles of the backward products,
 *   backwardlohi stores the second highest doubles of the backward products,
 *   backwardhilo stores the second lowest doubles of the backward products,
 *   backwardlolo stores the lowest doubles of the backward products,
 *                backward[nvr-3] contains the derivative with respect
 *                to the first variable idx[0] if nvr > 2;
 *   crosshihi    stores the highest doubles of the cross products,
 *   crosslohi    stores the second highest doubles of the cross products,
 *   crosshilo    stores the second lowest doubles of the cross products,
 *   crosslolo    stores the lowest doubles of the cross products,
 *                cross[k] contains the derivatve with respect to
 *                variable idx[k+1]. */

void GPU_cmplx4_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double *inputrehihi, double *inputrelohi,
   double *inputrehilo, double *inputrelolo,
   double *inputimhihi, double *inputimlohi,
   double *inputimhilo, double *inputimlolo,
   double *forwardrehihi, double *forwardrelohi,
   double *forwardrehilo, double *forwardrelolo,
   double *forwardimhihi, double *forwardimlohi,
   double *forwardimhilo, double *forwardimlolo,
   double *backwardrehihi, double *backwardrelohi,
   double *backwardrehilo, double *backwardrelolo,
   double *backwardimhihi, double *backwardimlohi,
   double *backwardimhilo, double *backwardimlolo,
   double *crossrehihi, double *crossrelohi,
   double *crossrehilo, double *crossrelolo,
   double *crossimhihi, double *crossimlohi,
   double *crossimhilo, double *crossimlolo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in quad double precision.
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
 *   cffrehihi      highest doubles of real parts of coefficients;
 *   cffrelohi      second highest doubles of real parts of coefficients;
 *   cffrehilo      second lowest doubles of real parts of coefficients;
 *   cffrelolo      lowest doubles of the real parts of the coefficients;
 *   cffimhihi      highest doubles of imaginary parts of the coefficients;
 *   cffimlohi      second highest doubles of imaginary parts of coefficients;
 *   cffimhilo      second lowest doubles of imaginary parts of coefficients;
 *   cffimlolo      lowest doubles of imaginary parts of coefficients;
 *   inputrehihi    holds the highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelohi    holds the second highest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrehilo    holds the second lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputrelolo    holds the lowest doubles of the real parts
 *                  of the input series for all variables;
 *   inputimhihi    holds the highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlohi    holds the second highest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimhilo    holds the second lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   inputimlolo    holds the lowest doubles of the imaginary parts
 *                  of the input series for all variables;
 *   forwardrehihi  holds work space for the highest doubles of nvr forward
 *                  products, for all real parts of the coefficients,
 *                  forwardrehi[k] holds space for deg+1 doubles;
 *   forwardrelohi  holds work space for the second highest doubles of nvr
 *                  forward products, for all real parts of the coefficients,
 *                  forwardrelohi[k] holds space for deg+1 doubles;
 *   forwardrehilo  holds work space for the second lowest doubles of nvr
 *                  forward products, for all real parts of the coefficients,
 *                  forwardrehilo[k] holds space for deg+1 doubles;
 *   forwardrelolo  holds work space for the lowest doubles of nvr forward
 *                  products, for all real parts of the coefficients,
 *                  forwardrelolo[k] holds space for deg+1 doubles;
 *   forwardimhihi  holds work space for the highest doubles of nvr forward
 *                  products, for all imaginary parts of the coefficients,
 *                  forwardimhihi[k] holds space for deg+1 doubles;
 *   forwardimlohi  holds work space for the second highest doubles of nvr
 *                  forward products, for imaginary parts of coefficients,
 *                  forwardimlohi[k] holds space for deg+1 doubles;
 *   forwardimhilo  holds work space for the second lowest doubles of nvr
 *                  forward products, for imaginary parts of coefficients,
 *                  forwardimhilo[k] holds space for deg+1 doubles;
 *   forwardimlolo  holds work space for the lowest doubles of nvr forward
 *                  products, for all imaginary parts of the coefficients,
 *                  forwardimlolo[k] holds space for deg+1 doubles;
 *   backwardrehihi holds work space for all highest doubles of nvr-2 backward
 *                  products, for all real parts of the coefficients,
 *                  backwardrehihi[k] holds space for deg+1 doubles;
 *   backwardrelohi holds work space for all second highest doubles of nvr-2
 *                  backward products, for all real parts of the coefficients,
 *                  backwardrelohi[k] holds space for deg+1 doubles;
 *   backwardrehilo holds work space for all second lowest doubles of nvr-2
 *                  backward products, for all real parts of the coefficients,
 *                  backwardrehilo[k] holds space for deg+1 doubles;
 *   backwardrelolo holds work space for all lowest doubles of nvr-2 backward
 *                  products, for all real parts of the coefficients,
 *                  backwardrelolo[k] holds space for deg+1 doubles;
 *   backwardimhihi holds work space for the highest doubles of nvr-2 backward
 *                  products, for all imaginary parts of the coefficients,
 *                  backwardimhihi[k] holds space for deg+1 doubles;
 *   backwardimlohi holds work space for the second highest doubles of nvr-2
 *                  backward products, for imaginary parts of coefficients,
 *                  backwardimlohi[k] holds space for deg+1 doubles;
 *   backwardimhilo holds work space for the second lowest doubles of nvr-2
 *                  backward products, for imaginary parts of coefficients,
 *                  backwardimhilo[k] holds space for deg+1 doubles;
 *   backwardimlolo holds work space for the lowest doubles of nvr-2 backward
 *                  products, for all imaginary parts of the coefficients,
 *                  backwardimlolo[k] holds space for deg+1 doubles;
 *   crossrehihi    holds work space for the highest doubles of nvr-2 cross
 *                  products,for the real parts of the coefficients,
 *                  crossrehihi[k] holds space for deg+1 doubles;
 *   crossrelohi    holds work space for the second highest doubles of nvr-2
 *                  cross products,for the real parts of the coefficients,
 *                  crossrelohi[k] holds space for deg+1 doubles;
 *   crossrehilo    holds work space for the second lowest doubles of nvr-2
 *                  cross products,for the real parts of the coefficients,
 *                  crossrehilo[k] holds space for deg+1 doubles;
 *   crossrelolo    holds work space for the lowest doubles of nvr-2 cross
 *                  products, for the real parts of the coefficients,
 *                  crossrelolo[k] holds space for deg+1 doubles;
 *   crossimhihi    holds work space for the highest doubles of nvr-2 cross
 *                  products, for the imaginary parts of the coefficients,
 *                  crossimhihi[k] holds space for deg+1 doubles.
 *   crossimlohi    holds work space for the second highest doubles of nvr-2
 *                  cross products, for the imaginary parts of coefficients,
 *                  crossimlohi[k] holds space for deg+1 doubles.
 *   crossimhilo    holds work space for the second lowest doubles of nvr-2
 *                  cross products, for the imaginary parts of coefficients,
 *                  crossimhilo[k] holds space for deg+1 doubles;
 *   crossimlolo    holds work space for the lowest doubles of nvr-2 cross
 *                  products, for the imaginary parts of the coefficients,
 *                  crossimlolo[k] holds space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrehihi  stores the highest doubles of the real parts
 *                  of the forward products,
 *   forwardrelohi  stores the second highest doubles of the real parts
 *                  of the forward products,
 *   forwardrehilo  stores the second lowest doubles of the real parts
 *                  of the forward products,
 *   forwardrelolo  stores the lowest doubles of the real parts
 *                  of the forward products,
 *   forwardimhihi  stores the highest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimlohi  stores the second highest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimhilo  stores the lowest doubles of the imaginary parts
 *                  of the forward products,
 *   forwardimlolo  stores the lowest doubles of the imaginary parts
 *                  of the forward products,
 *                  forward[nvr-1] contains the value of the product,
 *                  forward[nvr-2] contains the derivative with respect
 *                  to the last variable idx[nvr-1];
 *   backwardrehihi stores the highest doubles of the real parts
 *                  of the backward products,
 *   backwardrelohi stores the second highest doubles of the real parts
 *                  of the backward products,
 *   backwardrehilo stores the second lowest doubles of the real parts
 *                  of the backward products,
 *   backwardrelolo stores the lowest doubles of the real parts
 *                  of the backward products,
 *   backwardrehihi stores the highest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardrelohi stores the second highest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardrehilo stores the second lowest doubles of the imaginary parts
 *                  of the backward products,
 *   backwardrelolo stores the lowest doubles of the imaginary parts
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
 *                  cross[k] contains the derivatve with respect to
 *                  variable idx[k+1]. */

void GPU_dbl4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
 *   BS           number of threads in one block, must be deg+1;
 *   dim          total number of variables;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffhihi      deg+1 highest doubles of the coefficient series;
 *   cfflohi      deg+1 second highest doubles of the coefficient series;
 *   cffhilo      deg+1 second lowest doubles of the coefficient series;
 *   cfflolo      deg+1 lowest doubles of the coefficient series;
 *   inputhihi    holds the highest doubles of the input series
 *                for all variables in the monomial;
 *   inputlohi    holds the second highest doubles of the input series
 *                for all variables in the monomial;
 *   inputhilo    holds the second lowest doubles of the input series
 *                for all variables in the monomial;
 *   inputlolo    holds the lowest doubles of the input series
 *                for all variables in the monomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputhihi   stores the highest doubles of the derivatives and 
 *                the value of the monomial,
 *   outputlohi   stores the second highest doubles of the derivatives and 
 *                the value of the monomial,
 *   outputhilo   stores the second lowest doubles of the derivatives and 
 *                the value of the monomial,
 *   outputlolo   stores the lowest doubles of the derivatives and 
 *                the value of the monomial,
 *                output[idx[k]], for k from 0 to nvr, contains the
 *                deriviative with respect to the variable idx[k];
 *                output[dim] contains the value of the product. */

void GPU_cmplx4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
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
 *   BS           number of threads in one block, must be deg+1;
 *   dim          total number of variables;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffrehihi    highest doubles of the real parts of the coefficients;
 *   cffrelohi    second highest doubles of the real parts of the coefficients;
 *   cffrehilo    second lowest doubles of the real parts of the coefficients;
 *   cffrelolo    lowest doubles of the real parts of the coefficients;
 *   cffimhihi    highest doubles of the imag parts of the coefficients;
 *   cffimlohi    second highest doubles of the imag parts of the coefficients;
 *   cffimhilo    second lowest doubles of the imag parts of the coefficients;
 *   cffimlolo    lowest doubles of the imag parts of the coefficients;
 *   inputrehihi  holds the highest doubles of the real parts
 *                of the input series for all variables;
 *   inputrelohi  holds the second highest doubles of the real parts
 *                of the input series for all variables;
 *   inputrehilo  holds the second lowest doubles of the real parts
 *                of the input series for all variables;
 *   inputrelolo  holds the lowest doubles of the real parts
 *                of the input series for all variables;
 *   inputimhihi  holds the highest doubles of the imaginary parts
 *                of the input series for all variables;
 *   inputimlohi  holds the second highest doubles of the imaginary parts
 *                of the input series for all variables;
 *   inputimhilo  holds the second lowest doubles of the imaginary parts
 *                of the input series for all variables;
 *   inputimlolo  holds the lowest doubles of the imaginary parts
 *                of the input series for all variables;
 *   outputrehihi has space allocated for dim+1 series of degree deg, for the
 *                highest doubles of the real parts of the output;
 *   outputrelohi has space allocated for dim+1 series of degree deg, for the
 *                second highest doubles of the real parts of the output;
 *   outputrehilo has space allocated for dim+1 series of degree deg, for the
 *                second lowest doubles of the real parts of the output;
 *   outputrelolo has space allocated for dim+1 series of degree deg, for the
 *                lowest doubles of the real parts of the output;
 *   outputimhihi has space allocated for dim+1 series of degree deg, for the
 *                highest doubles of the imaginary parts of the output;
 *   outputimlohi has space allocated for dim+1 series of degree deg, for the
 *                second highest doubles of the imaginary parts of the output;
 *   outputimhilo has space allocated for dim+1 series of degree deg, for the
 *                second lowest doubles of the imaginary parts of the output;
 *   outputimlolo has space allocated for dim+1 series of degree deg, for the
 *                lowest doubles of the imaginary parts of the output.
 *
 * ON RETURN :
 *   outputrehihi stores the highest doubles of the real parts
 *                of the derivatives and the value,
 *   outputrelohi stores the second highest doubles of the real parts
 *                of the derivatives and the value,
 *   outputrehilo stores the second lowest doubles of the real parts
 *                of the derivatives and the value,
 *   outputrelolo stores the lowest doubles of the real parts
 *                of the derivatives and the value,
 *   outputimhihi stores the highest doubles of the imaginary parts
 *                of the derivatives and the value,
 *   outputimlohi stores the second highest doubles of the imaginary parts
 *                of the derivatives and the value,
 *   outputimhilo stores the lowest doubles of the imaginary parts
 *                of the derivatives and the value,
 *   outputimlolo stores the lowest doubles of the imaginary parts
 *                of the derivatives and the value,
 *                output[idx[k]], for k from 0 to nvr, contains the
 *                derivative with respect to the variable idx[k];
 *                output[dim] contains the value of the product. */

#endif
