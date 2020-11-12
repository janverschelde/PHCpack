/* The file dbl3_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in triple double precision. */

#ifndef __dbl3_monomials_host_h__
#define __dbl3_monomials_host_h__

void CPU_dbl3_speel
 ( int nvr, int deg, int *idx, double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **forwardhi, double **forwardmi, double **forwardlo,
   double **backwardhi, double **backwardmi, double **backwardlo,
   double **crosshi, double **crossmi, double **crosslo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   multiplied with a coefficient series of the same degree,
 *   for real coefficients in triple double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffhi      deg+1 doubles for the high doubles in the coefficient series;
 *   cffmi      deg+1 doubles for the middle doubles in the coefficient series;
 *   cfflo      deg+1 doubles for the low doubles in the coefficient series;
 *   inputhi    contains the high doubles of the coefficients of the series
 *              for all variables in the monomial;
 *   inputmi    contains the middle doubles of the coefficients of the series
 *              for all variables in the monomial;
 *   inputlo    contains the low doubles of the coefficients of the series
 *              for all variables in the monomial;
 *   forwardhi  contains work space for the high doubles of nvr forward 
 *              products, forwardhi[k] contains space for deg+1 doubles;
 *   forwardmi  contains work space for the middle doubles of nvr forward 
 *              products, forwardmi[k] contains space for deg+1 doubles;
 *   forwardlo  contains work space for the low doubles of nvr forward 
 *              products, forwardlo[k] contains space for deg+1 doubles;
 *   backwardhi contains work space for the high doubles of nvr-2 backward
 *              products, backwardhi[k] contains space for deg+1 doubles;
 *   backwardmi contains work space for the middle doubles of nvr-2 backward
 *              products, backwardmi[k] contains space for deg+1 doubles;
 *   backwardlo contains work space for the low doubles of nvr-2 backward
 *              products, backwardlo[k] contains space for deg+1 doubles;
 *   crosshi    contains work space for the high doubles of nvr-2 cross
 *              products, crosshi[k] contains space for deg+1 doubles;
 *   crossmi    contains work space for the middle doubles of nvr-2 cross
 *              products, crosshi[k] contains space for deg+1 doubles;
 *   crosslo    contains work space for the low doubles of nvr-2 cross
 *              products, crosslo[k] contains space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardhi  accumulates the high doubles of the forward products,
 *   forwardmi  accumulates the middle doubles of the forward products,
 *   forwardlo  accumulates the low doubles of the forward products,
 *              forward[nvr-1] contains the value of the product,
 *              forward[nvr-2] contains the derivative with respect
 *              to the last variable idx[nvr-1];
 *   backwardhi accumulates the high doubles of the backward products,
 *   backwardmi accumulates the middle doubles of the backward products,
 *   backwardlo accumulates the low doubles of the backward products,
 *              backward[nvr-3] contains the derivative with respect
 *              to the first variable idx[0];
 *   crosshi    stores the high doubles of the cross products,
 *   crossmi    stores the middle doubles of the cross products,
 *   crosslo    stores the low doubles of the cross products,
 *              cross[k] contains the derivatve with respect to
 *              variable idx[k+1]. */

void CPU_cmplx3_speel
 ( int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **forwardrehi, double **forwardremi, double **forwardrelo,
   double **forwardimhi, double **forwardimmi, double **forwardimlo,
   double **backwardrehi, double **backwardremi, double **backwardrelo,
   double **backwardimhi, double **backwardimmi, double **backwardimlo,
   double **crossrehi, double **crossremi, double **crossrelo,
   double **crossimhi, double **crossimmi, double **crossimlo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in triple double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr          number of variables in the product;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffrehi      high doubles of the real parts of the coefficients
 *                of the series of the product;
 *   cffremi      middle doubles of the real parts of the coefficients
 *                of the series of the product;
 *   cffrelo      low doubles of the real parts of the coefficients
 *                of the series of the product;
 *   cffimhi      high doubles of the imaginary pars of the coefficient
 *                of the series of the product;
 *   cffimmi      middle doubles of the imaginary pars of the coefficient
 *                of the series of the product;
 *   cffimlo      low doubles of the imaginary pars of the coefficient
 *                of the series of the product;
 *   inputrehi    contains the high doubles of the real parts of the
 *                coefficients of the series for all variables;
 *   inputremi    contains the middle doubles of the real parts of the
 *                coefficients of the series for all variables;
 *   inputrelo    contains the low doubles of the real parts of the
 *                coefficients of the series for all variables;
 *   inputimhi    contains the high doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   inputimmi    contains the middle doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   inputimlo    contains the low doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   forwardrehi  contains work space for the high doubles of nvr forward
 *                products, forwardrehi[k] can hold deg+1 doubles;
 *   forwardremi  contains work space for the middle doubles of nvr forward
 *                products, forwardremi[k] can hold deg+1 doubles;
 *   forwardrelo  contains work space for the low doubles of nvr forward
 *                products, forwardrelo[k] can hold deg+1 doubles;
 *   forwardimhi  contains work space for the high doubles of nvr forward
 *                products, forwardimhi[k] can hold deg+1 doubles;
 *   forwardimmi  contains work space for the middle doubles of nvr forward
 *                products, forwardimmi[k] can hold deg+1 doubles;
 *   forwardimlo  contains work space for the low doubles of nvr forward
 *                products, forwardimlo[k] can hold deg+1 doubles;
 *   backwardrehi contains work space for the high doubles of nvr-2 backward
 *                products, backwardrehi[k] can hold deg+1 doubles;
 *   backwardremi contains work space for the middle doubles of nvr-2 backward
 *                products, backwardremi[k] can hold deg+1 doubles;
 *   backwardrelo contains work space for the low doubles of nvr-2 backward
 *                products, backwardrelo[k] can hold deg+1 doubles;
 *   backwardimhi contains work space for the high doubles of nvr-2 backward
 *                products, backwardimhi[k] can hold deg+1 doubles;
 *   backwardimmi contains work space for the middle doubles of nvr-2 backward
 *                products, backwardimmi[k] can hold deg+1 doubles;
 *   backwardimlo contains work space for the low doubles of nvr-2 backward
 *                products, backwardimlo[k] can hold deg+1 doubles;
 *   crossrehi    contains work space for the high doubles of nvr-2 cross
 *                products, crossrehi[k] can hold deg+1 doubles;
 *   crossremi    contains work space for the middle doubles of nvr-2 cross
 *                products, crossremi[k] can hold deg+1 doubles;
 *   crossrelo    contains work space for the low doubles of nvr-2 cross
 *                products, crossrelo[k] can hold deg+1 doubles;
 *   crossimhi    contains work space for the high doubles of nvr-2 cross
 *                products, crossimhi[k] can hold deg+1 doubles;
 *   crossimmi    contains work space for the middle doubles of nvr-2 cross
 *                products, crossimmi[k] can hold deg+1 doubles;
 *   crossimlo    contains work space for the low doubles of nvr-2 cross
 *                products, crossimlo[k] can hold deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrehi  accumulates the high doubles of the real parts
 *                of the forward products,
 *   forwardremi  accumulates the middle doubles of the real parts
 *                of the forward products,
 *   forwardrelo  accumulates the low doubles of the real parts
 *                of the forward products,
 *   forwardimhi  accumulates the high doubles of the imaginary parts
 *                of the forward products,
 *   forwardimmi  accumulates the middle doubles of the imaginary parts
 *                of the forward products,
 *   forwardimlo  accumulates the low doubles of the imaginary parts
 *                of the forward products,
 *                forward[nvr-1] contains the value of the product,
 *                forward[nvr-2] contains the derivative with respect
 *                to the last variable idx[nvr-1];
 *   backwardrehi accumulates the high doubles of the real parts
 *                of the backward products,
 *   backwardremi accumulates the middle doubles of the real parts
 *                of the backward products,
 *   backwardrelo accumulates the low doubles of the real parts
 *                of the backward products,
 *   backwardimhi accumulates the high doubles of the imaginary parts
 *                of the backward products,
 *   backwardimmi accumulates the middle doubles of the imaginary parts
 *                of the backward products,
 *   backwardimlo accumulates the low doubles of the imaginary parts
 *                of the backward products,
 *                backward[nvr-3] contains the derivative with respect
 *                to the first variable idx[0];
 *   crossrehi    stores the high doubles of the real parts
 *                of the cross products,
 *   crossremi    stores the middle doubles of the real parts
 *                of the cross products,
 *   crossrelo    stores the low doubles of the real parts
 *                of the cross products,
 *   crossimhi    stores the high doubles of the imaginary parts
 *                of the cross products,
 *   crossimmi    stores the middle doubles of the imaginary parts
 *                of the cross products,
 *   crossimlo    stores the low doubles of the imaginary parts
 *                of the cross products,
 *                cross[k] contains the derivatve with respect to
 *                the variable idx[k+1]. */

void CPU_dbl3_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   deg      truncation degree of the series;
 *   idx      as many indices as the value of nvr,
 *            idx[k] defines the place of the k-th variable,
 *            with input values in input[idx[k]];
 *   cffhi    high doubles of the coefficient series of the product;
 *   cffmi    middle doubles of the coefficient series of the product;
 *   cfflo    low doubles of the coefficient series of the product;
 *   inputhi  contains the high doubles of the coefficients of the series
 *            for all variables in the monomial;
 *   inputmi  contains the middle doubles of the coefficients of the series
 *            for all variables in the monomial;
 *   inputlo  contains the low doubles of the coefficients of the series
 *            for all variables in the monomial;
 *   outputhi has space allocated for dim+1 series of degree deg;
 *   outputmi has space allocated for dim+1 series of degree deg;
 *   outputlo has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputhi contains high doubles of the derivatives and the value,
 *   outputmi contains middle doubles of the derivatives and the value,
 *   outputlo contains low doubles of the derivatives and the value,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

void CPU_cmplx3_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffrehi    high doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffremi    middle doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffrelo    low doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffimhi    high doubles of the imaginary pars of the coefficients
 *              of the series of the product;
 *   cffimmi    middle doubles of the imaginary pars of the coefficients
 *              of the series of the product;
 *   cffimlo    low doubles of the imaginary pars of the coefficients
 *              of the series of the product;
 *   inputrehi  contains the high doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputremi  contains the middle doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputrelo  contains the low doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputimhi  contains the high doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimmi  contains the middle doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimlo  contains the low doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   outputrehi has space allocated for dim+1 series of degree deg;
 *   outputremi has space allocated for dim+1 series of degree deg;
 *   outputrelo has space allocated for dim+1 series of degree deg;
 *   outputimhi has space allocated for dim+1 series of degree deg;
 *   outputimmi has space allocated for dim+1 series of degree deg;
 *   outputimlo has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputrehi contains the high doubles of the real parts
 *              of the derivatives and the value,
 *   outputremi contains the middle doubles of the real parts
 *              of the derivatives and the value,
 *   outputimlo contains the low doubles of the imaginary parts
 *              of the derivatives and the value,
 *              output[idx[k]], for k from 0 to nvr, contains the
 *              deriviative with respect to the variable idx[k];
 *              output[dim] contains the value of the product. */

#endif
