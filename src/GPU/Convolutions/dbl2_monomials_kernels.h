// The file dbl2_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in double double precision.

#ifndef __dbl2_monomials_kernels_h__
#define __dbl2_monomials_kernels_h__

void GPU_dbl2_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhi, double *cfflo, double *inputhi,
   double *inputlo, double *forwardhi, double *forwardlo, double *backwardhi,
   double *backwardlo, double *crosshi, double *crosslo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in double double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS         number of threads in one block, must be deg+1;
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffhi      deg+1 doubles for high doubles of the coefficient series;
 *   cfflo      deg+1 doubles for low doubles of the coefficient series;
 *   inputhi    stores the high doubles of the coefficients of the series
 *              for all variables in the monomial;
 *   inputlo    stores the high doubles of the coefficients of the series
 *              for all variables in the monomial;
 *   forwardhi  is work space for the high doubles of all nvr forward
 *              products, forwardhi has space for nvr*(deg+1) doubles;
 *   forwardlo  is work space for the low doubles of all nvr forward
 *              products, forwardlo has space for nvr*(deg+1) doubles;
 *   backwardhi is work space for the high doubles of all nvr-1 backward
 *              products, backwardhi has space for (nvr-1)*(deg+1) doubles;
 *   backwardlo is work space for the low doubles of all nvr-1 backward
 *              products, backwardlo has space for (nvr-1)*(deg+1) doubles;
 *   crosshi    is work space for the high doubles of all nvr-2 cross
 *              products, crosshi has space for (nvr-2)*(deg+1) doubles;
 *   crosslo    is work space for the low doubles of all nvr-2 cross
 *              products, crosslo has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardhi  stores the high doubles of the forward products,
 *   forwardlo  stores the low doubles of the forward products,
 *              forward[nvr-1] is the value of the product (as double**),
 *              forward[nvr-2] is the derivative with respect
 *              to the last variable idx[nvr-1] if nvr > 2;
 *   backwardhi stores the high doubles of the backward products,
 *   backwardlo stores the low doubles of the backward products,
 *              backward[nvr-2] is the derivative with respect
 *              to the first variable idx[0] if nvr > 2;
 *   crosshi    stores the high doubles of the cross products,
 *   crosslo    stores the low doubles of the cross products,
 *              cross[k] is the derivatve with respect to
 *              variable idx[k+1]. */

void GPU_cmplx2_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo,
   double *inputrehi, double *inputrelo, double *inputimhi, double *inputimlo,
   double *forwardrehi, double *forwardrelo, double *forwardimhi,
   double *forwardimlo, double *backwardrehi, double *backwardrelo,
   double *backwardimhi, double *backwardimlo, double *crossrehi,
   double *crossrelo, double *crossimhi, double *crossimlo );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in double double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cffre, cffim, inputre, and inputim are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forwardre, forwardim, backwardre, backwardim, crossre, and crossim
 *   are allocated on the device.
 *
 * ON ENTRY :
 *   BS           number of threads in one block, must be deg+1;
 *   nvr          number of variables in the product;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffrehi      high doubles of the real parts of the series coefficient
 *                the product;
 *   cffrelo      low doubles of the real parts of the series coefficient
 *                of the product;
 *   cffimhi      high doubles of the imaginary pars of the series coefficient
 *                of the product;
 *   cffimlo      low doubles of the imaginary pars of the series coefficient
 *                of the product;
 *   inputrehi    stores the high doubles of the real parts of the
 *                coefficients of the series for all variables;
 *   inputrelo    stores the low doubles of the real parts of the
 *                coefficients of the series for all variables;
 *   inputimhi    stores the high doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   inputimlo    stores the low doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   forwardrehi  is work space for the high doubles of nvr forward
 *                products, for all real parts of the coefficients,
 *                forwardrehi has space for nvr*(deg+1) doubles;
 *   forwardrelo  is work space for the low doubles of nvr forward
 *                products, for all real parts of the coefficients,
 *                forwardrelo has space for nvr*(deg+1) doubles;
 *   forwardimhi  is work space for the high doubles of nvr forward
 *                products, for all imaginary parts of the coefficients,
 *                forwardimhi has space for nvr*(deg+1) doubles;
 *   forwardimlo  is work space for the low doubles of nvr forward
 *                products, for all imaginary parts of the coefficients,
 *                forwardimlo has space for nvr*(deg+1) doubles;
 *   backwardrehi is work space for all high doubles of nvr-1 backward
 *                products, for all real parts of the coefficients,
 *                backwardrehi has space for (nvr-1)*(deg+1) doubles;
 *   backwardrelo is work space for all low doubles of nvr-1 backward
 *                products, for all real parts of the coefficients,
 *                backwardrelo has space for (nvr-1)*(deg+1) doubles;
 *   backwardimhi is work space for the high doubles of nvr-1 backward
 *                products, for all imaginary parts of the coefficients,
 *                backwardimhi has space for (nvr-1)*(deg+1) doubles;
 *   backwardimlo is work space for the low doubles of nvr-1 backward
 *                products, for all imaginary parts of the coefficients,
 *                backwardimlo has space for (nvr-1)*(deg+1) doubles;
 *   crossrehi    is work space for the high doubles of nvr-2 cross
 *                products,for the real parts of the coefficients,
 *                crossrehi has space for (nvr-2)*(deg+1) doubles;
 *   crossrelo    is work space for the low doubles of nvr-2 cross
 *                products,for the real parts of the coefficients,
 *                crossrelo has space for (nvr-2)*(deg+1) doubles;
 *   crossimhi    is work space for the high doubles of nvr-2 cross
 *                products, for the imaginary parts of the coefficients,
 *                crossimhi has space for (nvr-2)*(deg+1) doubles.
 *   crossimlo    is work space for the low doubles of nvr-2 cross
 *                products, for the imaginary parts of the coefficients,
 *                crossimlo has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardrehi  stores the high doubles of the real parts
 *                of the forward products,
 *   forwardrelo  stores the low doubles of the real parts
 *                of the forward products,
 *   forwardimhi  stores the high doubles of the imaginary parts
 *                of the forward products,
 *   forwardimlo  stores the low doubles of the imaginary parts
 *                of the forward products, as double**,
 *                forward[nvr-1] is the value of the product,
 *                forward[nvr-2] is the derivative with respect
 *                to the last variable idx[nvr-1];
 *   backwardrehi stores the high doubles of the real parts
 *                of the backward products,
 *   backwardrelo stores the low doubles of the real parts
 *                of the backward products,
 *   backwardimhi stores the high doubles of the imaginary parts of 
 *                the backward products,
 *   backwardimlo stores the low doubles of the imaginary parts of 
 *                the backward products,
 *                backward[nvr-2] is the derivative with respect
 *                to the first variable idx[0];
 *   crossrehi    stores the high doubles of the real parts
 *                of the cross products,
 *   crossrelo    stores the low doubles of the real parts
 *                of the cross products,
 *   crossimhi    stores the high doubles of the imaginary parts
 *                of the cross products,
 *   crossimlo    stores the low doubles of the imaginary parts
 *                of the cross products,
 *                cross[k] is the derivatve with respect to
 *                variable idx[k+1]. */

void GPU_dbl2_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cffhi, double *cfflo,
   double **inputhi, double **inputlo, double **outputhi, double **outputlo );
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
 *   cffhi    high doubles of the coefficient series of the product;
 *   cfflo    low doubles of the coefficient series of the product;
 *   inputhi  contains the high doubles of the coefficients of the series
 *            for all variables in the monomial;
 *   inputlo  contains the low doubles of the coefficients of the series
 *            for all variables in the monomial;
 *   outputhi has space allocated for dim+1 series of degree deg;
 *   outputlo has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputhi contains the high doubles of the derivatives and 
 *            the value of the monomial,
 *   outputlo contains the low doubles of the derivatives and 
 *            the value of the monomial,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

void GPU_cmplx2_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo,
   double **inputrehi, double **inputrelo, double **inputimhi,
   double **inputimlo, double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   BS         number of threads in one block, must be deg+1;
 *   dim        total number of variables;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffrehi    high doubles of the real parts of the coefficient series
 *              of the product;
 *   cffrelo    low doubles of the real parts of the coefficient series
 *              of the product;
 *   cffimhi    high doubles of the imaginary parts of the coefficient series
 *              of the product;
 *   cffimlo    low doubles of the imaginary parts of the coefficient series
 *              of the product;
 *   inputrehi  contains the high doubles of the real parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputrelo  contains the low doubles of the real parts of the coefficients
 *              of the series for all variables in the monomial;
 *   inputimhi  contains the high doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimlo  contains the low doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   outputrehi has space allocated for dim+1 series of degree deg, for the
 *              high doubles of the real parts of the output coefficients;
 *   outputrelo has space allocated for dim+1 series of degree deg, for the
 *              low doubles of the real parts of the output coefficients;
 *   outputimhi has space allocated for dim+1 series of degree deg, for the
 *              high doubles of the imaginary parts of the output;
 *   outputimlo has space allocated for dim+1 series of degree deg, for the
 *              low doubles of the imaginary parts of the output.
 *
 * ON RETURN :
 *   outputrehi contains the high doubles of the real parts of the
 *              derivatives and the value,
 *   outputrelo contains the low doubles of the real parts of the
 *              derivatives and the value,
 *   outputimhi contains the high doubles of the imaginary parts of the 
 *              derivatives and the value,
 *   outputimlo contains the low doubles of the imaginary parts of the 
 *              derivatives and the value,
 *              output[idx[k]], for k from 0 to nvr, contains the
 *              derivative with respect to the variable idx[k];
 *              output[dim] contains the value of the product. */

#endif
