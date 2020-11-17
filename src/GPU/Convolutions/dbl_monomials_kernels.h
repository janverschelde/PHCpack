// The file dbl_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in double precision.

#ifndef __dbl_monomials_kernels_h__
#define __dbl_monomials_kernels_h__

void GPU_dbl_speel
 ( int BS, int nvr, int deg, int *idx, double *cff, double *input,
   double *forward, double *backward, double *cross );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS       number of threads in one block, must be deg+1;
 *   nvr      number of variables in the product;
 *   deg      truncation degree of the series;
 *   idx      as many indices as the value of nvr,
 *            idx[k] defines the place of the k-th variable,
 *            with input values in input[idx[k]];
 *   cff      deg+1 doubles for the coefficient series of the monomial;
 *   input    stores the coefficients of the power series
 *            for all variables in the monomial,
 *   forward  is work space for all nvr forward products,
 *            forward has space for nvr*(deg+1) doubles;
 *   backward is work space for all nvr-1 backward products;
 *            backward has space for (nvr-1)*(deg+1) doubles;
 *   cross    is work space for all nvr-2 cross products;
 *            cross has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forward  stores the forward products,
 *            forward[nvr-1] is the value of the product (as double**),
 *            forward[nvr-2] is the derivative with respect
 *            to the last variable idx[nvr-1];
 *   backward stores the backward products,
 *            backward[nvr-2] contains the derivative with respect
 *            to the first variable idx[0];
 *   cross    stores the cross products,
 *            cross[k] contains the derivative with respect to
 *            variable idx[k+1]. */

void GPU_cmplx_speel
 ( int BS, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double *inputre, double *inputim, double *forwardre,
   double *forwardim, double *backwardre, double *backwardim,
   double *crossre, double *crossim );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cffre, cffim, inputre, and inputim are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forwardre, forwardim, backwardre, backwardim, crossre, and crossim
 *   are allocated on the device.
 *
 * ON ENTRY :
 *   BS         number of threads in one block, must be deg+1;
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffre      real parts of the series coefficient of the product;
 *   cffim      imaginary pars of the series coefficients of the product;
 *   inputre    stores the real parts of the coefficients of the series
 *              for all variables in the monomial;
 *   inputim    stores the imaginary parts of the coefficients
 *              of the series for all variables in the monomial;
 *   forwardre  is work space for all nvr forward products,
 *              for all real parts of the coefficients,
 *              forwardre has space for nvr*(deg+1) doubles;
 *   forwardim  is work space for all nvr forward products,
 *              for all imaginary parts of the coefficients,
 *              forwardim has space for nvr*(deg+1) doubles;
 *   backwardre is work space for all nvr-1 backward products,
 *              for all real parts of the coefficients,
 *              backwardre has space for (nvr-1)*(deg+1) doubles;
 *   backwardim is work space for all nvr-1 backward products,
 *              for all imaginary parts of the coefficients,
 *              backwardim has space for (nvr-1)*(deg+1) doubles;
 *   crossre    is work space for all nvr-2 cross products,
 *              for the real parts of the coefficients,
 *              crossre has space for (nvr-2)*(deg+1) doubles;
 *   crossim    is work space for all nvr-2 cross products,
 *              for the imaginary parts of the coefficients,
 *              crossim has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardre  stores the real parts of the forward products,
 *   forwardim  stores the imaginary parts of the forward products,
 *              forward[nvr-1] is the value of the product (as double**)
 *              forward[nvr-2] is the derivative with respect
 *              to the last variable idx[nvr-1];
 *   backwardre stores the real parts of the backward products,
 *   backwardre stores the imaginary parts of the backward products,
 *              backward[nvr-2] is the derivative with respect
 *              to the first variable idx[0];
 *   crossre    stores the real parts of the cross products,
 *   crossim    stores the imaginary parts of the cross products,
 *              cross[k] contains the derivatve with respect to
 *              variable idx[k+1]. */

void GPU_dbl_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cff,
   double **input, double **output );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   BS       number of threads in one block, must be deg+1;
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

void GPU_cmplx_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **outputre, double **outputim );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   BS       number of threads in one block, must be deg+1;
 *   dim      total number of variables;
 *   deg      truncation degree of the series;
 *   idx      as many indices as the value of nvr,
 *            idx[k] defines the place of the k-th variable,
 *            with input values in input[idx[k]];
 *   cffre    real parts of the coefficient series of the product;
 *   cffim    imaginary parts of the coefficient series of the product;
 *   inputre  contains the real parts of the coefficients of the series
 *            for all variables in the monomial;
 *   inputim  contains the imaginary parts of the coefficients of the 
 *            series for all variables in the monomial;
 *   outputre has space allocated for dim+1 series of degree deg,
 *            for the real parts of the output coefficients;
 *   outputim has space allocated for dim+1 series of degree deg,
 *            for the imaginary parts of the output coefficients.
 *
 * ON RETURN :
 *   outputre contains real parts of the derivatives and the value,
 *   outputim contains imaginary parts of the derivatives and the value,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

#endif
