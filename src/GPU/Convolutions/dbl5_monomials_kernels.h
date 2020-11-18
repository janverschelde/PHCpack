// The file dbl5_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in penta double precision.

#ifndef __dbl5_monomials_kernels_h__
#define __dbl5_monomials_kernels_h__

void GPU_dbl5_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double *inputtb, double *inputix, double *inputmi,
   double *inputrg, double *inputpk,
   double *forwardtb, double *forwardix, double *forwardmi,
   double *forwardrg, double *forwardpk,
   double *backwardtb, double *backwardix, double *backwardmi,
   double *backwardrg, double *backwardpk,
   double *crosstb, double *crossix, double *crossmi,
   double *crossrg, double *crosspk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in penta double precision.
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
 *   cfftb      deg+1 highest doubles of the coefficient series;
 *   cffix      deg+1 second highest doubles of the coefficient series;
 *   cffmi      deg+1 middle doubles of the coefficient series;
 *   cffrg      deg+1 second lowest doubles of the coefficient series;
 *   cffpk      deg+1 lowest doubles of the coefficient series;
 *   inputtb    stores the highest doubles of the input series
 *              for all variables in the monomial;
 *   inputix    stores the second highest doubles of the input series
 *              for all variables in the monomial;
 *   inputmi    stores the middle doubles of the input series
 *              for all variables in the monomial;
 *   inputrg    stores the second lowest doubles of the input series
 *              for all variables in the monomial;
 *   inputpk    stores the lowest doubles of the input series
 *              for all variables in the monomial;
 *   forwardtb  is work space for the highest doubles of nvr 
 *              forward products, has space for nvr*(deg+1) doubles;
 *   forwardix  is work space for the second highest doubles of nvr 
 *              forward products, has space for nvr*(deg+1) doubles;
 *   forwardmi  is work space for the middle doubles of nvr
 *              forward products, has space for nvr*(deg+1) doubles;
 *   forwardrg  is work space for the second lowest doubles of nvr
 *              forward products, has space for nvr*(deg+1) doubles;
 *   forwardpk  is work space for the lowest doubles of nvr
 *              forward products, has space for nvr*(deg+1) doubles;
 *   backwardtb is work space for the highest doubles of nvr-2
 *              backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardix is work space for the second highest doubles of nvr-2
 *              backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardmi is work space for the middle doubles of nvr-2
 *              backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardrg is work space for the second lowest doubles of nvr-2
 *              backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardpk is work space for the lowest doubles of nvr-2
 *              backward products, has space for (nvr-2)*(deg+1) doubles;
 *   crosstb    is work space for the highest doubles of nvr-2
 *              cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossix    is work space for the second highest doubles of nvr-2
 *              cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossmi    is work space for the middle doubles of nvr-2
 *              cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossrg    is work space for the second lowest doubles of nvr-2
 *              cross products, has space for (nvr-2)*(deg+1) doubles.
 *   crosspk    is work space for the lowest doubles of nvr-2
 *              cross products, has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardtb  stores the highest doubles of the forward products,
 *   forwardix  stores the second highest doubles of the forward products,
 *   forwardmi  stores the middle doubles of the forward products,
 *   forwardrg  stores the second lowest doubles of the forward products,
 *   forwardpk  stores the lowest doubles of the forward products,
 *              forward[nvr-1] contains the value of the product,
 *              forward[nvr-2] contains the derivative with respect
 *              to the last variable idx[nvr-1] if nvr > 2;
 *   backwardtb stores the highest doubles of the backward products,
 *   backwardix stores the second highest doubles of the backward products,
 *   backwardmi stores the middle doubles of the backward products,
 *   backwardrg stores the second lowest doubles of the backward products,
 *   backwardpk stores the lowest doubles of the backward products,
 *              backward[nvr-3] contains the derivative with respect
 *              to the first variable idx[0] if nvr > 2;
 *   crosstb    stores the highest doubles of the cross products,
 *   crossix    stores the second highest doubles of the cross products,
 *   crossmi    stores the middle doubles of the cross products,
 *   crossrg    stores the second lowest doubles of the cross products,
 *   crosspk    stores the lowest doubles of the cross products,
 *              cross[k] contains the derivative with respect to
 *              variable idx[k+1]. */

void GPU_cmplx5_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double *inputretb, double *inputreix, double *inputremi,
   double *inputrerg, double *inputrepk,
   double *inputimtb, double *inputimix, double *inputimmi,
   double *inputimrg, double *inputimpk,
   double *forwardretb, double *forwardreix, double *forwardremi,
   double *forwardrerg, double *forwardrepk,
   double *forwardimtb, double *forwardimix, double *forwardimmi,
   double *forwardimrg, double *forwardimpk,
   double *backwardretb, double *backwardreix, double *backwardremi,
   double *backwardrerg, double *backwardrepk,
   double *backwardimtb, double *backwardimix, double *backwardimmi,
   double *backwardimrg, double *backwardimpk,
   double *crossretb, double *crossreix, double *crossremi,
   double *crossrerg, double *crossrepk,
   double *crossimtb, double *crossimix, double *crossimmi,
   double *crossimrg, double *crossimpk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in penta double precision.
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
 *   cffretb      highest doubles of the real parts
 *                of the series coefficient of the product;
 *   cffreix      second highest doubles of the real parts
 *                of the series coefficient of the product;
 *   cffremi      middle doubles of the real parts
 *                of the series coefficient of the product;
 *   cffrerg      second lowest doubles of the real parts
 *                of the series coefficient of the product;
 *   cffrepk      lowest doubles of the real parts
 *                of the series coefficient of the product;
 *   cffimtb      highest doubles of the imaginary parts
 *                of the series coefficient of the product;
 *   cffimix      second highest doubles of the imaginary parts
 *                of the series coefficient of the product;
 *   cffimmi      middle doubles of the imaginary parts
 *                of the series coefficient of the product;
 *   cffimrg      second lowest doubles of the imaginary parts
 *                of the series coefficient of the product;
 *   cffimpk      lowest doubles of the imaginary parts
 *                of the series coefficient of the product;
 *   inputretb    stores the highest doubles of the real parts
 *                of the coefficients of the input series;
 *   inputreix    stores the second highest doubles of the real parts
 *                of the coefficients of the input series;
 *   inputremi    stores the middle doubles of the real parts
 *                of the coefficients of the input series;
 *   inputrerg    stores the second lowest doubles of the real parts
 *                of the coefficients of the input series;
 *   inputrepk    stores the lowest doubles of the real parts
 *                of the coefficients of the input series;
 *   inputimtb    stores the highest doubles of the imaginary parts
 *                of the coefficients of the input series;
 *   inputimix    stores the second highest doubles of the imaginary parts
 *                of the coefficients of the input series;
 *   inputimmi    stores the middle doubles of the imaginary parts
 *                of the coefficients of the input series;
 *   inputimrg    stores the second lowest doubles of the imaginary parts
 *                of the coefficients of the input series;
 *   inputimpk    stores the lowest doubles of the imaginary parts
 *                of the coefficients of the input series;
 *   forwardretb  is work space for the highest doubles
 *                for all real parts of the nvr forward products,
 *                forwardretb has space for nvr*(deg+1) doubles;
 *   forwardreix  is work space for the second highest doubles
 *                for all real parts of the nvr forward products,
 *                forwardreix has space for nvr*(deg+1) doubles;
 *   forwardremi  is work space for the middle doubles
 *                for all real parts of the nvr forward products,
 *                forwardremi has space for nvr*(deg+1) doubles;
 *   forwardrerg  is work space for the second lowest doubles
 *                for all real parts of the nvr forward products,
 *                forwardrerg has space for nvr*(deg+1) doubles;
 *   forwardrepk  is work space for the lowest doubles
 *                for all real parts of the nvr forward products,
 *                forwardrepk has space for nvr*(deg+1) doubles;
 *   forwardimtb  is work space for the highest doubles
 *                for all imaginary parts of the nvr forward products,
 *                forwardimtb has space for nvr*(deg+1) doubles;
 *   forwardimix  is work space for the second highest doubles
 *                for all imaginary parts of the nvr forward products,
 *                forwardimix has space for nvr*(deg+1) doubles;
 *   forwardimmi  is work space for the middle doubles
 *                for all imaginary parts of the nvr forward products,
 *                forwardimmi has space for nvr*(deg+1) doubles;
 *   forwardimrg  is work space for the second lowest doubles
 *                for all imaginary parts of the nvr forward products,
 *                forwardimrg is space for nvr*(deg+1) doubles;
 *   forwardimpk  is work space for the lowest doubles
 *                for all imaginary parts of the nvr forward products,
 *                forwardimpk is space for nvr*(deg+1) doubles;
 *   backwardretb is work space for all highest doubles
 *                for all real parts of the nvr-2 backward products,
 *                backwardretb has space for (nvr-2)*(deg+1) doubles;
 *   backwardreix is work space for all second highest doubles
 *                for all real parts of the nvr-2 backward products,
 *                backwardreix has space for (nvr-2)*(deg+1) doubles;
 *   backwardremi is work space for all middle doubles
 *                for all real parts of the nvr-2 backward products,
 *                backwardremi has space for (nvr-2)*(deg+1) doubles;
 *   backwardrerg is work space for all second lowest doubles
 *                for all real parts of the nvr-2 backward products,
 *                backwardrerg has space for (nvr-2)*(deg+1) doubles;
 *   backwardrepk is work space for all lowest doubles
 *                for all real parts of the nvr-2 backward products,
 *                backwardrepk has space for (nvr-2)*(deg+1) doubles;
 *   backwardimtb is work space for the highest doubles
 *                for all imaginary parts of the nvr-2 backward products,
 *                backwardimtb has space for (nvr-2)*(deg+1) doubles;
 *   backwardimix is work space for the second highest doubles
 *                for all imaginary parts of the nvr-2 backward products,
 *                backwardimix has space for (nvr-2)*(deg+1) doubles;
 *   backwardimmi is work space for the middle doubles
 *                for all imaginary parts of the nvr-2 backward products,
 *                backwardimmi has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrg is work space for the second lowest doubles
 *                for all imaginary parts of the nvr-2 backward products,
 *                backwardimrg has space for (nvr-2)*(deg+1) doubles;
 *   backwardimpk is work space for the lowest doubles
 *                for all imaginary parts of the nvr-2 backward products,
 *                backwardimpk has space for (nvr-2)*(deg+1) doubles;
 *   crossretb    is work space for the highest doubles
 *                for the real parts of the nvr-2 cross products,
 *                crossretb has space for (nvr-2)*(deg+1) doubles;
 *   crossreix    is work space for the second highest doubles
 *                for the real parts of the nvr-2 cross products,
 *                crossreix has space for (nvr-2)*(deg+1) doubles;
 *   crossremi    is work space for the middle doubles
 *                for the real parts of the nvr-2 cross products,
 *                crossremi has space for (nvr-2)*(deg+1) doubles;
 *   crossrerg    is work space for the second lowest doubles
 *                for the real parts of nvr-2 cross products,
 *                crossrerg has space for (nvr-2)*(deg+1) doubles;
 *   crossrepk    is work space for the lowest doubles
 *                for the real parts of the nvr-2 cross products,
 *                crossrepk has space for (nvr-2)*(deg+1) doubles;
 *   crossimtb    is work space for the highest doubles
 *                for the imaginary parts of nvr-2 cross products,
 *                crossimtb has space for (nvr-2)*(deg+1) doubles;
 *   crossimix    is work space for the second highest doubles
 *                for the imaginary parts of nvr-2 cross products,
 *                crossimix has space for (nvr-2)*(deg+1) doubles;
 *   crossimmi    is work space for the middle doubles
 *                for the imaginary parts of nvr-2 cross products,
 *                crossimmi has space for (nvr-2)*(deg+1) doubles;
 *   crossimrg    is work space for the second lowest doubles
 *                for the imaginary parts of nvr-2 cross products,
 *                crossimrg has space for (nvr-2)*(deg+1) doubles;
 *   crossimpk    is work space for the lowest doubles
 *                for the imaginary parts of nvr-2 cross products,
 *                crossimpk has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardretb  stores the highest doubles of the real parts
 *                of the forward products,
 *   forwardreix  stores the second highest doubles of the real parts
 *                of the forward products,
 *   forwardremi  stores the middle doubles of the real parts
 *                of the forward products,
 *   forwardrerg  stores the second lowest doubles of the real parts
 *                of the forward products,
 *   forwardrepk  stores the lowest doubles of the real parts
 *                of the forward products,
 *   forwardimtb  stores the highest doubles of the imaginary parts
 *                of the forward products,
 *   forwardimix  stores the second highest doubles of the imaginary parts
 *                of the forward products,
 *   forwardimmi  stores the middle doubles of the imaginary parts
 *                of the forward products,
 *   forwardimrg  stores the second lowest doubles of the imaginary parts
 *                of the forward products,
 *   forwardimpk  stores the lowest doubles of the imaginary parts
 *                forward[nvr-1] contains the value of the product,
 *                forward[nvr-2] contains the derivative with respect
 *                to the last variable idx[nvr-1];
 *   backwardretb stores the highest doubles of the real parts
 *                of the backward products,
 *   backwardreix stores the second highest doubles of the real parts
 *                of the backward products,
 *   backwardremi stores the middle doubles of the real parts
 *                of the backward products,
 *   backwardrerg stores the second lowest doubles of the real parts
 *                of the backward products,
 *   backwardrepk stores the lowest doubles of the real parts
 *                of the backward products,
 *   backwardretb stores the highest doubles of the imaginary parts
 *                of the backward products,
 *   backwardreix stores the second highest doubles of the imaginary parts
 *                of the backward products,
 *   backwardremi stores the middle doubles of the imaginary parts
 *                of the backward products,
 *   backwardrerg stores the second lowest doubles of the imaginary parts
 *                of the backward products,
 *   backwardrepk stores the lowest doubles of the imaginary parts of 
 *                the backward products,
 *                backward[nvr-3] contains the derivative with respect
 *                to the first variable idx[0];
 *   crossretb    stores the highest doubles of the real parts
 *                of the cross products,
 *   crossreix    stores the second highest doubles of the real parts
 *                of the cross products,
 *   crossremi    stores the middle doubles of the real parts
 *                of the cross products,
 *   crossrerg    stores the second lowest doubles of the real parts
 *                of the cross products,
 *   crossrepk    stores the lowest doubles of the real parts
 *                of the cross products,
 *   crossimtb    stores the highest doubles of the imaginary parts
 *                of the cross products,
 *   crossimix    stores the second highest doubles of the imaginary parts
 *                of the cross products,
 *   crossimmi    stores the middle doubles of the imaginary parts
 *                of the cross products,
 *   crossimrg    stores the second lowest doubles of the imaginary parts
 *                of the cross products,
 *   crossimpk    stores the lowest doubles of the imaginary parts
 *                of the cross products,
 *                cross[k] contains the derivative with respect to
 *                variable idx[k+1]. */

void GPU_dbl5_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk );
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
 *   cfftb    deg+1 highest doubles of the coefficient series;
 *   cffix    deg+1 second highest doubles of the coefficient series;
 *   cffmi    deg+1 middle doubles of the coefficient series;
 *   cffrg    deg+1 second lowest doubles of the coefficient series;
 *   cffpk    deg+1 lowest doubles of the coefficient series;
 *   inputtb  stores the highest doubles of the input series
 *            for all variables in the monomial;
 *   inputix  stores the second highest doubles of the input series
 *            for all variables in the monomial;
 *   inputmi  stores the middle doubles of the input series
 *            for all variables in the monomial;
 *   inputrg  stores the second lowest doubles of the input series
 *            for all variables in the monomial;
 *   inputpk  stores the lowest doubles of the input series
 *            for all variables in the monomial;
 *   outputtb has space allocated for dim+1 series of degree deg;
 *   outputix has space allocated for dim+1 series of degree deg;
 *   outputmi has space allocated for dim+1 series of degree deg;
 *   outputrg has space allocated for dim+1 series of degree deg;
 *   outputpk has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputtb stores the highest doubles of the derivatives and 
 *            the value of the monomial,
 *   outputix stores the second highest doubles of the derivatives and 
 *            the value of the monomial,
 *   outputmi stores the middle doubles of the derivatives and 
 *            the value of the monomial,
 *   outputrg stores the second lowest doubles of the derivatives and 
 *            the value of the monomial,
 *   outputpk stores the lowest doubles of the derivatives and 
 *            the value of the monomial,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

void GPU_cmplx5_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk );
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
 *   cffretb    highest doubles of the real parts
 *              of the series coefficient of the product;
 *   cffreix    second highest doubles of the real parts
 *              of the series coefficient of the product;
 *   cffremi    middle doubles of the real parts
 *              of the series coefficient of the product;
 *   cffrerg    second lowest doubles of the real parts
 *              of the series coefficient of the product;
 *   cffrepk    lowest doubles of the real parts
 *              of the series coefficient of the product;
 *   cffimtb    highest doubles of the imaginary parts
 *              of the series coefficient of the product;
 *   cffimix    second highest doubles of the imaginary parts
 *              of the series coefficient of the product;
 *   cffimmi    middle doubles of the imaginary parts
 *              of the series coefficient of the product;
 *   cffimrg    second lowest doubles of the imaginary parts
 *              of the series coefficient of the product;
 *   cffimpk    lowest doubles of the imaginary parts
 *              of the series coefficient of the product;
 *   inputretb  stores the highest doubles of the real parts
 *              of the coefficients of the input series;
 *   inputreix  stores the second highest doubles of the real parts
 *              of the coefficients of the input series;
 *   inputremi  stores the middle doubles of the real parts
 *              of the coefficients of the input series;
 *   inputrerg  stores the second lowest doubles of the real parts
 *              of the coefficients of the input series;
 *   inputrepk  stores the lowest doubles of the real parts
 *              of the coefficients of the input series;
 *   inputimtb  stores the highest doubles of the imaginary parts
 *              of the coefficients of the input series;
 *   inputimix  stores the second highest doubles of the imaginary parts
 *              of the coefficients of the input series;
 *   inputimmi  stores the middle doubles of the imaginary parts
 *              of the coefficients of the input series;
 *   inputimrg  stores the second lowest doubles of the imaginary parts
 *              of the coefficients of the input series;
 *   inputimpk  stores the lowest doubles of the imaginary parts
 *              of the coefficients of the input series;
 *   outputretb has space allocated for dim+1 series of degree deg,
 *              for the real highest doubles of the output;
 *   outputreix has space allocated for dim+1 series of degree deg,
 *              for the real second highest doubles of the output;
 *   outputremi has space allocated for dim+1 series of degree deg,
 *              for the real middle doubles of the output;
 *   outputrerg has space allocated for dim+1 series of degree deg,
 *              for the real second lowest doubles of the output;
 *   outputrepk has space allocated for dim+1 series of degree deg,
 *              for the real lowest doubles of the output;
 *   outputimtb has space allocated for dim+1 series of degree deg,
 *              for the imaginary highest doubles of the output;
 *   outputimix has space allocated for dim+1 series of degree deg,
 *              for the imaginary second highest doubles of the output;
 *   outputimmi has space allocated for dim+1 series of degree deg,
 *              for the imaginary middle doubles of the output;
 *   outputimrg has space allocated for dim+1 series of degree deg,
 *              for the imaginary second lowest doubles of the output;
 *   outputimpk has space allocated for dim+1 series of degree deg,
 *              for the imaginary lowest doubles of the output.
 *
 * ON RETURN :
 *   outputretb stores the highest doubles of the real parts
 *              of the derivatives and the value,
 *   outputreix stores the second highest doubles of the real parts
 *              of the derivatives and the value,
 *   outputremi stores the middle doubles of the real parts
 *              of the derivatives and the value,
 *   outputrerg stores the second lowest doubles of the real parts
 *              of the derivatives and the value,
 *   outputrepk stores the lowest doubles of the real parts
 *              of the derivatives and the value,
 *   outputimtb stores the highest doubles of the imaginary parts
 *              of the derivatives and the value,
 *   outputimix stores the second highest doubles of the imaginary parts
 *              of the derivatives and the value,
 *   outputimmi stores the middle doubles of the imaginary parts
 *              of the derivatives and the value,
 *   outputimrg stores the second lowest doubles of the imaginary parts
 *              derivatives and the value,
 *   outputimpk stores the lowest doubles of the imaginary parts
 *              derivatives and the value,
 *              output[idx[k]], for k from 0 to nvr, contains the
 *              derivative with respect to the variable idx[k];
 *              output[dim] contains the value of the product. */

#endif
