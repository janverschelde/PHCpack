/* The file dbl_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in double precision. */

#ifndef __dbl_monomials_host_h__
#define __dbl_monomials_host_h__

void CPU_dbl_speel
 ( int nvr, int deg, int *idx, double *cff, double **input,
   double **forward, double **backward, double **cross,
   bool verbose=false, int monidx=0 );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   multiplied with a coefficient series of the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr      number of variables in the product;
 *   deg      truncation degree of the series;
 *   idx      as many indices as the value of nvr,
 *            idx[k] defines the place of the k-th variable,
 *            with input values in input[idx[k]];
 *   cff      deg+1 doubles for the coefficient series of the monomial;
 *   input    contains the coefficients of the power series
 *            for all variables in the monomial;
 *   forward  contains work space for all nvr forward products,
 *            forward[k] contains space for deg+1 doubles;
 *   backward contains work space for all nvr-2 backward products;
 *            backward[k] contains space for deg+1 doubles;
 *   cross    contains work space for all nvr-2 cross products;
 *            cross[k] contains space for deg+1 doubles;
 *   verbose  if true, writes one line to screen for every convolution;
 *   monidx   index of the monomial, only needed if verbose.
 *
 * ON RETURN :
 *   forward  accumulates the forward products,
 *            forward[nvr-1] contains the value of the product,
 *            forward[nvr-2] contains the derivative with respect
 *            to the last variable idx[nvr-1];
 *   backward accumulates the backward products,
 *            backward[nvr-3] contains the derivative with respect
 *            to the first variable idx[0];
 *   cross    stores the cross products,
 *            cross[k] contains the derivatve with respect to
 *            variable idx[k+1]. */

void CPU_cmplx_speel
 ( int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **forwardre,
   double **forwardim, double **backwardre, double **backwardim,
   double **crossre, double **crossim );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffre      real parts of the coefficients of the series of the product;
 *   cffim      imaginary pars of the coefficient of the series of the product;
 *   inputre    contains the real parts of the coefficients of the series
 *              for all variables in the monomial;
 *   inputim    contains the imaginary parts of the coefficients of the series
 *              for all variables in the monomial;
 *   forwardre  contains work space for all nvr forward products,
 *              forwardre[k] can hold deg+1 doubles for real parts;
 *   forwardim  contains work space for all nvr forward products,
 *              forwardim[k] can hold deg+1 doubles for imaginary parts;
 *   backwardre contains work space for all nvr-2 backward products;
 *              backwardre[k] can hold deg+1 doubles for real parts;
 *   backwardim contains work space for all nvr-2 backward products;
 *              backwardim[k] can hold deg+1 doubles for imaginary parts;
 *   crossre    contains work space for all nvr-2 cross products;
 *              crossre[k] can hold deg+1 doubles for real parts;
 *   crossim    contains work space for all nvr-2 cross products;
 *              crossim[k] can hold deg+1 doubles for imaginary parts.
 *
 * ON RETURN :
 *   forwardre  accumulates the real parts of the forward products,
 *   forwardim  accumulates the imaginary parts of the forward products,
 *              forward[nvr-1] contains the value of the product,
 *              forward[nvr-2] contains the derivative with respect
 *              to the last variable idx[nvr-1];
 *   backwardre accumulates the real parts of the backward products,
 *   backwardim accumulates the imaginary parts of the backward products,
 *              backward[nvr-3] contains the derivative with respect
 *              to the first variable idx[0];
 *   crossre    stores the real parts of the cross products,
 *   crossim    stores the imaginary parts of the cross products,
 *              cross[k] contains the derivatve with respect to
 *              the variable idx[k+1]. */

void CPU_dbl_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cff,
   double **input, double **output );
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
 *   cff      coefficient power series of the product;
 *   input    contains the coefficients of the power series
 *            for all variables in the monomial;
 *   output   space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the monomial,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

void CPU_cmplx_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **outputre, double **outputim );
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
 *   cffre    real parts of the coefficients of the series of the product;
 *   cffim    imaginary pars of the coefficients of the series of the product;
 *   inputre  contains the real parts of the coefficients of the series
 *            for all variables in the monomial;
 *   inputim  contains the imaginary parts of the coefficients of the series
 *            for all variables in the monomial;
 *   outputre has space allocated for dim+1 series of degree deg;
 *   outputim has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputre contains real parts of the derivatives and the value,
 *   outputim contains imaginary parts of the derivatives and the value,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

#endif
