/* The file dbl5_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in penta double precision. */

#ifndef __dbl5_monomials_host_h__
#define __dbl5_monomials_host_h__

void CPU_dbl5_speel
 ( int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi, double **inputrg,
   double **inputpk, double **forwardtb, double **forwardix,
   double **forwardmi, double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk, double **crosstb,
   double **crossix, double **crossmi, double **crossrg, double **crosspk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   multiplied with a coefficient series of the same degree,
 *   for real coefficients in penta double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr        number of variables in the product;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cfftb      deg+1 highest doubles in the coefficient series;
 *   cffix      deg+1 second highest doubles in the coefficient series;
 *   cffmi      deg+1 middle doubles in the coefficient series;
 *   cffrg      deg+1 second lowest doubles in the coefficient series;
 *   cffpk      deg+1 lowest doubles in the coefficient series;
 *   inputhi    holds the highest doubles of the input series
 *              for all variables in the monomial;
 *   inputix    holds the second highest doubles of the input series
 *              for all variables in the monomial;
 *   inputmi    holds the middle doubles of the input series
 *              for all variables in the monomial;
 *   inputrg    holds the second lowest doubles of the input series
 *              for all variables in the monomial;
 *   inputpk    holds the lowest doubles of the input series
 *              for all variables in the monomial;
 *   forwardtb  is work space for the highest doubles of nvr forward 
 *              products, forwardtb[k] has space for deg+1 doubles;
 *   forwardix  is work space for the second highest doubles of nvr forward 
 *              products, forwardix[k] has space for deg+1 doubles;
 *   forwardmi  is work space for the middle doubles of nvr forward 
 *              products, forwardmi[k] has space for deg+1 doubles;
 *   forwardrg  is work space for the second lowest doubles of nvr forward 
 *              products, forwardrg[k] has space for deg+1 doubles;
 *   forwardpk  is work space for the lowest doubles of nvr forward 
 *              products, forwardpk[k] has space for deg+1 doubles;
 *   backwardtb is work space for the highest doubles of nvr-2 backward
 *              products, backwardtb[k] has space for deg+1 doubles;
 *   backwardix is work space for the second highest doubles of nvr-2 backward
 *              products, backwardix[k] has space for deg+1 doubles;
 *   backwardmi is work space for the middle doubles of nvr-2 backward
 *              products, backwardmi[k] has space for deg+1 doubles;
 *   backwardrg is work space for the second lowest doubles of nvr-2 backward
 *              products, backwardrg[k] has space for deg+1 doubles;
 *   backwardpk is work space for the lowest doubles of nvr-2 backward
 *              products, backwardpk[k] has space for deg+1 doubles;
 *   crosstb    is work space for the highest doubles of nvr-2 cross
 *              products, crosstb[k] has space for deg+1 doubles;
 *   crossix    is work space for the second highest doubles of nvr-2 cross
 *              products, crossix[k] has space for deg+1 doubles;
 *   crossmi    is work space for the middle doubles of nvr-2 cross
 *              products, crossmi[k] has space for deg+1 doubles;
 *   crossrg    is work space for the second lowest doubles of nvr-2 cross
 *              products, crossrg[k] has space for deg+1 doubles;
 *   crosspk    is work space for the lowest doubles of nvr-2 cross
 *              products, crosspk[k] has space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardtb  stores the highest doubles of the forward products,
 *   forwardix  stores the second highest doubles of the forward products,
 *   forwardmi  stores the middle doubles of the forward products,
 *   forwardrg  stores the second lowest doubles of the forward products,
 *   forwardpk  stores the lowest doubles of the forward products,
 *              forward[nvr-1] contains the value of the product,
 *              forward[nvr-2] contains the derivative with respect
 *              to the last variable idx[nvr-1];
 *   backwardtb stores the highest doubles of the backward products,
 *   backwardix stores the second highest doubles of the backward products,
 *   backwardmi stores the middle doubles of the backward products,
 *   backwardrg stores the second lowest doubles of the backward products,
 *   backwardpk stores the lowest doubles of the backward products,
 *              backward[nvr-3] contains the derivative with respect
 *              to the first variable idx[0];
 *   crosstb    stores the highest doubles of the cross products,
 *   crossix    stores the second highest doubles of the cross products,
 *   crossmi    stores the middle doubles of the cross products,
 *   crossrg    stores the second lowest doubles of the cross products,
 *   crosspk    stores the lowest doubles of the cross products,
 *              cross[k] contains the derivatve with respect to
 *              variable idx[k+1]. */

void CPU_cmplx5_speel
 ( int nvr, int deg, int *idx, double *cffretb, double *cffreix,
   double *cffremi, double *cffrerg, double *cffrepk, double *cffimtb,
   double *cffimix, double *cffimmi, double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **forwardretb, double **forwardreix, double **forwardremi,
   double **forwardrerg, double **forwardrepk,
   double **forwardimtb, double **forwardimix, double **forwardimmi,
   double **forwardimrg, double **forwardimpk,
   double **backwardretb, double **backwardreix, double **backwardremi,
   double **backwardrerg, double **backwardrepk,
   double **backwardimtb, double **backwardimix, double **backwardimmi,
   double **backwardimrg, double **backwardimpk,
   double **crossretb, double **crossreix, double **crossremi,
   double **crossrerg, double **crossrepk,
   double **crossimtb, double **crossimix, double **crossimmi,
   double **crossimrg, double **crossimpk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in penta double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr          number of variables in the product;
 *   deg          truncation degree of the series;
 *   idx          as many indices as the value of nvr,
 *                idx[k] defines the place of the k-th variable,
 *                with input values in input[idx[k]];
 *   cffretb      highest doubles of the real parts of the coefficients
 *                of the series of the product;
 *   cffreix      second highest doubles of the real parts of the
 *                coefficients of the series of the product;
 *   cffremi      middle doubles of the real parts of the coefficients
 *                of the series of the product;
 *   cffrerg      second lowest doubles of the real parts of the
 *                coefficients of the series of the product;
 *   cffrepk      lowest doubles of the real parts of the
 *                coefficients of the series of the product;
 *   cffimtb      highest doubles of the imaginary parts of the coefficients
 *                of the series of the product;
 *   cffimix      second highest doubles of the imaginary parts of the
 *                coefficients of the series of the product;
 *   cffimmi      middle doubles of the imaginary parts of the coefficients
 *                of the series of the product;
 *   cffimrg      second lowest doubles of the imaginary parts of the
 *                coefficients of the series of the product;
 *   cffimpk      lowest doubles of the imaginary parts of the coefficients
 *                of the series of the product;
 *   inputretb    holds the highest doubles of the real parts of the 
 *                coefficients of the series for all variables;
 *   inputreix    holds the second highest doubles of the real parts of the 
 *                coefficients of the series for all variables;
 *   inputremi    holds the middle doubles of the real parts of the 
 *                coefficients of the series for all variables;
 *   inputrerg    holds the second lowest doubles of the real parts of the 
 *                coefficients of the series for all variables;
 *   inputrepk    holds the lowest doubles of the real parts of the 
 *                coefficients of the series for all variables;
 *   inputimtb    holds the highest doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   inputimix    holds the second highest doubles of the imaginary parts of
 *                the coefficients of the series for all variables;
 *   inputimmi    holds the middle doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   inputimrg    holds the second lowest doubles of the imaginary parts of
 *                the coefficients of the series for all variables;
 *   inputimpk    holds the lowest doubles of the imaginary parts of the
 *                coefficients of the series for all variables;
 *   forwardretb  is work space for the highest doubles of nvr forward
 *                products, forwardretb[k] has space for deg+1 doubles;
 *   forwardreix  is work space for the second highest doubles of nvr forward
 *                products, forwardreix[k] has space for deg+1 doubles;
 *   forwardremi  is work space for the middle doubles of nvr forward
 *                products, forwardremi[k] has space for deg+1 doubles;
 *   forwardrerg  is work space for the second lowest doubles of nvr forward
 *                products, forwardrerg[k] has space for deg+1 doubles;
 *   forwardrepk  is work space for the lowest doubles of nvr forward
 *                products, forwardrepk[k] has space for deg+1 doubles;
 *   forwardimtb  is work space for the highest doubles of nvr forward
 *                products, forwardimtb[k] has space for deg+1 doubles;
 *   forwardimix  is work space for the second highest doubles of nvr forward
 *                products, forwardimix[k] has space for deg+1 doubles;
 *   forwardimmi  is work space for the middle doubles of nvr forward
 *                products, forwardimmi[k] has space for deg+1 doubles;
 *   forwardimrg  is work space for the second lowest doubles of nvr forward
 *                products, forwardimrg[k] has space for deg+1 doubles;
 *   forwardimpk  is work space for the lowest doubles of nvr forward
 *                products, forwardimpk[k] has space for deg+1 doubles;
 *   backwardretb is work space for the highest doubles of nvr-2 backward
 *                products, backwardretb[k] has space for deg+1 doubles;
 *   backwardreix is work space for second highest doubles of nvr-2 backward
 *                products, backwardreix[k] has space for deg+1 doubles;
 *   backwardremi is work space for the middle doubles of nvr-2 backward
 *                products, backwardremi[k] has space for deg+1 doubles;
 *   backwardrerg is work space for second lowest doubles of nvr-2 backward
 *                products, backwardrerg[k] has space for deg+1 doubles;
 *   backwardrepk is work space for the lowest doubles of nvr-2 backward
 *                products, backwardrepk[k] has space for deg+1 doubles;
 *   backwardimtb is work space for the highest doubles of nvr-2 backward
 *                products, backwardimtb[k] has space for deg+1 doubles;
 *   backwardimix is work space for second highest doubles of nvr-2 backward
 *                products, backwardimix[k] has space for deg+1 doubles;
 *   backwardimmi is work space for the middle doubles of nvr-2 backward
 *                products, backwardimmi[k] has space for deg+1 doubles;
 *   backwardimrg is work space for second lowest doubles of nvr-2 backward
 *                products, backwardimrg[k] has space for deg+1 doubles;
 *   backwardimpk is work space for the lowest doubles of nvr-2 backward
 *                products, backwardimpk[k] has space for deg+1 doubles;
 *   crossretb    is work space for the highest doubles of nvr-2 cross
 *                products, crossretb[k] has space for deg+1 doubles;
 *   crossreix    is work space for the second highest doubles of nvr-2 cross
 *                products, crossreix[k] has space for deg+1 doubles;
 *   crossremi    is work space for the middle doubles of nvr-2 cross
 *                products, crossremi[k] has space for deg+1 doubles;
 *   crossrerg    is work space for the second lowest doubles of nvr-2 cross
 *                products, crossrerg[k] has space for deg+1 doubles;
 *   crossrepk    is work space for the lowest doubles of nvr-2 cross
 *                products, crossrepk[k] has space for deg+1 doubles;
 *   crossimtb    is work space for the highest doubles of nvr-2 cross
 *                products, crossimtb[k] has space for deg+1 doubles;
 *   crossimix    is work space for the second highest doubles of nvr-2 cross
 *                products, crossimix[k] has space for deg+1 doubles;
 *   crossimmi    is work space for the middle doubles of nvr-2 cross
 *                products, crossimmi[k] has space for deg+1 doubles;
 *   crossimrg    is work space for the second lowest doubles of nvr-2 cross
 *                products, crossimrg[k] has space for deg+1 doubles.
 *   crossimpk    is work space for the lowest doubles of nvr-2 cross
 *                products, crossimpk[k] has space for deg+1 doubles.
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
 *                of the forward products,
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
 *   backwardimtb stores the highest doubles of the imaginary parts
 *                of the backward products,
 *   backwardimix stores the second highest doubles of the imaginary parts
 *                of the backward products,
 *   backwardimmi stores the middle doubles of the imaginary parts
 *                of the backward products,
 *   backwardimrg stores the second lowest doubles of the imaginary parts
 *                of the backward products,
 *   backwardimpk stores the lowest doubles of the imaginary parts
 *                of the backward products,
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
 *                cross[k] contains the derivatve with respect to
 *                the variable idx[k+1]. */

void CPU_dbl5_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi, double **inputrg,
   double **inputpk, double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk );
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
 *   cfftb    deg+1 highest doubles in the coefficient series;
 *   cffix    deg+1 second highest doubles in the coefficient series;
 *   cffmi    deg+1 middle doubles in the coefficient series;
 *   cffrg    deg+1 second lowest doubles in the coefficient series;
 *   cffpk    deg+1 lowest doubles in the coefficient series;
 *   inputtb  holds the highest doubles of the input series
 *            for all variables in the monomial;
 *   inputix  holds the second highest doubles of the input series
 *            for all variables in the monomial;
 *   inputmi  holds the middle doubles of the input series
 *            for all variables in the monomial;
 *   inputrg  holds the second lowest doubles of the input series
 *            for all variables in the monomial;
 *   inputpk  holds the lowest doubles of the input series
 *            for all variables in the monomial;
 *   outputtb has space allocated for dim+1 series of degree deg;
 *   outputix has space allocated for dim+1 series of degree deg;
 *   outputmi has space allocated for dim+1 series of degree deg;
 *   outputrg has space allocated for dim+1 series of degree deg;
 *   outputpk has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputtb stores highest doubles of the derivatives and the value,
 *   outputix stores second highest doubles of the derivatives and the value,
 *   outputmi stores middle doubles of the derivatives and the value,
 *   outputrg stores second lowest doubles of the derivatives and the value,
 *   outputpk stores lowest doubles of the derivatives and the value,
 *            output[idx[k]], for k from 0 to nvr, contains the
 *            deriviative with respect to the variable idx[k];
 *            output[dim] contains the value of the product. */

void CPU_cmplx5_evaldiff
 ( int dim, int nvr, int deg, int *idx, double *cffretb, double *cffreix,
   double *cffremi, double *cffrerg, double *cffrepk, double *cffimtb,
   double *cffimix, double *cffimmi, double *cffimrg, double *cffimpk,
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
 *   dim        total number of variables;
 *   deg        truncation degree of the series;
 *   idx        as many indices as the value of nvr,
 *              idx[k] defines the place of the k-th variable,
 *              with input values in input[idx[k]];
 *   cffretb    highest doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffreix    second highest doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffremi    middle doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffrerg    second lowest doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffrepk    lowest doubles of the real parts of the coefficients
 *              of the series of the product;
 *   cffimtb    highest doubles of the imaginary parts of the coefficients
 *              of the series of the product;
 *   cffimix    second highest doubles of the imaginary parts of the
 *              coefficients of the series of the product;
 *   cffimmi    middle doubles of the imaginary parts of the coefficients
 *              of the series of the product;
 *   cffimrg    second lowest doubles of the imaginary parts of the
 *              coefficients of the series of the product;
 *   cffimpk    lowest doubles of the imaginary parts of the coefficients
 *              of the series of the product;
 *   inputretb  holds the highest doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputreix  holds the second highest doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputremi  holds the middle doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputrerg  holds the second lowest doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputrepk  holds the lowest doubles of the real parts of the 
 *              coefficients of the series for all variables in the monomial;
 *   inputimtb  holds the highest doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimix  holds the second highest doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimmi  holds the middle doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimpk  holds the second lowest doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   inputimpk  holds the lowest doubles of the imaginary parts of the
 *              coefficients of the series for all variables in the monomial;
 *   outputretb has space allocated for dim+1 series of degree deg;
 *   outputreix has space allocated for dim+1 series of degree deg;
 *   outputremi has space allocated for dim+1 series of degree deg;
 *   outputrerg has space allocated for dim+1 series of degree deg;
 *   outputrepk has space allocated for dim+1 series of degree deg;
 *   outputimtb has space allocated for dim+1 series of degree deg;
 *   outputimix has space allocated for dim+1 series of degree deg;
 *   outputimmi has space allocated for dim+1 series of degree deg;
 *   outputimrg has space allocated for dim+1 series of degree deg;
 *   outputimpk has space allocated for dim+1 series of degree deg.
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
 *              of the derivatives and the value,
 *   outputimpk stores the lowest doubles of the imaginary parts
 *              of the derivatives and the value,
 *              output[idx[k]], for k from 0 to nvr, contains the
 *              deriviative with respect to the variable idx[k];
 *              output[dim] contains the value of the product. */

#endif
