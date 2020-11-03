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

__global__ void GPU_dbl_speel
 ( int nvr, int deg, int *idx, double *input,
   double *forward, double *backward, double *cross );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED : nvr >= 2.
 *
 * ON ENTRY :
 *   nvr      number of variables in the product;
 *   deg      truncation degree of the series;
 *   idx      as many indices as the value of nvr,
 *            idx[k] defines the place of the k-th variable,
 *            with input values in input[idx[k]];
 *   input    contains the coefficients of the power series
 *            for all variables in the monomial;
 *   forward  contains work space for all nvr-1 forward products,
 *            forward[k] contains space for deg doubles;
 *   backward contains work space for all nvr-2 backward products;
 *            backward[k] contains space for deg doubles;
 *   cross    contains work space for all nvr-2 cross products;
 *            cross[k] contains space for deg doubles.
 *
 * ON RETURN :
 *   forward  accumulates the forward products,
 *            forward[nvr-2] contains the value of the product,
 *            forward[nvr-3] contains the derivative with respect
 *            to the last variable idx[nvr-1] if nvr > 2;
 *   backward accumulates the backward products,
 *            backward[0] contains the derivative with respect
 *            to the first variable idx[0] if nvr > 2;
 *   cross    stores the cross products,
 *            cross[k] contains the derivatve with respect to
 *            variable idx[k+1]. */

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

#endif
