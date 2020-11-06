// The file dbl_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in double precision.

#ifndef __dbl_monomials_kernels_h__
#define __dbl_monomials_kernels_h__

__device__ void dbl_convolute
 ( double *x, double *y, double *z, int dim, int k );
/*
 * DESCRIPTION :
 *   Thread k returns in z[k] the k-th component of the convolution
 *   of x and y.  All vectors are of dimension dim. */

__device__ void cmplx_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim, int k );
/*
 * DESCRIPTION :
 *   Thread k returns in z[k] the k-th component of the convolution
 *   of x and y, with real and imaginary parts in re and im arrays.
 *   All arrays are of dimension dim. */

__global__ void GPU_dbl_speel
 ( int nvr, int deg, int *idx, double *cff, double *input,
   double *forward, double *backward, double *cross );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
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
 *            cross[k] contains space for deg+1 doubles.
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

__global__ void GPU_cmplx_speel
 ( int nvr, int deg, int *idx, double *cffre, double *cffim,
   double *inputre, double *inputim, double *forwardre,
   double *forwardim, double *backwardre, double *backwardim,
   double *crossre, double *crossim );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffre      real parts of the series coefficient of the product;
 *   cffim      imaginary pars of the series coefficients of the product;
 *   inputre    contains the real parts of the coefficients of the series
 *              for all variables in the monomial;
 *   inputim    contains the imaginary parts of the coefficients
 *              of the series for all variables in the monomial;
 *   forwardre  contains work space for all nvr forward products,
 *              for all real parts of the coefficients,
 *              forwardre[k] contains space for deg+1 doubles;
 *   forwardim  contains work space for all nvr forward products,
 *              for all imaginary parts of the coefficients,
 *              forwardim[k] contains space for deg+1 doubles;
 *   backwardre contains work space for all nvr-2 backward products,
 *              for all real parts of the coefficients,
 *              backwardre[k] contains space for deg+1 doubles;
 *   backwardim contains work space for all nvr-2 backward products,
 *              for all imaginary parts of the coefficients,
 *              backwardim[k] contains space for deg+1 doubles;
 *   crossre    contains work space for all nvr-2 cross products,
 *              for the real parts of the coefficients,
 *              crossre[k] contains space for deg+1 doubles;
 *   crossim    contains work space for all nvr-2 cross products,
 *              for the imaginary parts of the coefficients,
 *              crossim[k] contains space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardre  accumulates the real parts of the forward products,
 *   forwardim  accumulates the imaginary parts of the forward products,
 *              forward[nvr-1] contains the value of the product,
 *              forward[nvr-2] contains the derivative with respect
 *              to the last variable idx[nvr-1];
 *   backwardre accumulates the real parts of the backward products,
 *   backwardre accumulates the imaginary parts of the backward products,
 *              backward[nvr-3] contains the derivative with respect
 *              to the first variable idx[0]
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
