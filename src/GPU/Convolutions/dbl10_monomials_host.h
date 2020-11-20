/* The file dbl10_monomials_host.h specifies functions to evaluate and
 * differentiate a monomial at power series truncated to the same degree,
 * in deca double precision. */

#ifndef __dbl10_monomials_host_h__
#define __dbl10_monomials_host_h__

void CPU_dbl10_speel
 ( int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **forwardrtb, double **forwardrix, double **forwardrmi,
   double **forwardrrg, double **forwardrpk,
   double **forwardltb, double **forwardlix, double **forwardlmi,
   double **forwardlrg, double **forwardlpk,
   double **backwardrtb, double **backwardrix, double **backwardrmi,
   double **backwardrrg, double **backwardrpk,
   double **backwardltb, double **backwardlix, double **backwardlmi,
   double **backwardlrg, double **backwardlpk,
   double **crossrtb, double **crossrix, double **crossrmi,
   double **crossrrg, double **crossrpk,
   double **crossltb, double **crosslix, double **crosslmi,
   double **crosslrg, double **crosslpk );
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
 *   nvr         number of variables in the product;
 *   deg         truncation degree of the series;
 *   idx         as many indices as the value of nvr,
 *               idx[k] defines the place of the k-th variable,
 *               with input values in input[idx[k]];
 *   cffrtb      deg+1 higest doubles in the coefficient series;
 *   cffrix      deg+1 second higest doubles in the coefficient series;
 *   cffrmi      deg+1 third highest doubles in the coefficient series;
 *   cffrrg      deg+1 fourth highest doubles in the coefficient series;
 *   cffrpk      deg+1 fifth highest doubles in the coefficient series;
 *   cffltb      deg+1 fifth lowest doubles in the coefficient series;
 *   cfflix      deg+1 fourth lowest doubles in the coefficient series;
 *   cfflmi      deg+1 third lowest doubles in the coefficient series;
 *   cfflrg      deg+1 second lowest doubles in the coefficient series;
 *   cfflpk      deg+1 lowest doubles in the coefficient series;
 *   inputrtb    holds the highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrix    holds the second highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrmi    holds the third highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrrg    holds the fourth highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrpk    holds the fifth highest doubles of the input series
 *               for all variables in the monomial;
 *   inputltb    holds the fifth lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlix    holds the fourth lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlmi    holds the third lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlrg    holds the second lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlpk    holds the lowest doubles of the input series
 *               for all variables in the monomial;
 *   forwardrtb  is work space for the highest doubles of nvr forward 
 *               products, forwardrtb[k] has space for deg+1 doubles;
 *   forwardrix  is work space for the second highest doubles of nvr forward 
 *               products, forwardrix[k] has space for deg+1 doubles;
 *   forwardrmi  is work space for the third highest doubles of nvr forward 
 *               products, forwardrmi[k] has space for deg+1 doubles;
 *   forwardrrg  is work space for the fourth highest doubles of nvr forward 
 *               products, forwardrrg[k] has space for deg+1 doubles;
 *   forwardrpk  is work space for the fifth highest doubles of nvr forward 
 *               products, forwardrpk[k] has space for deg+1 doubles;
 *   forwardltb  is work space for the fifth lowest doubles of nvr forward 
 *               products, forwardltb[k] has space for deg+1 doubles;
 *   forwardlix  is work space for the fourth lowest doubles of nvr forward 
 *               products, forwardlix[k] has space for deg+1 doubles;
 *   forwardlmi  is work space for the third lowest doubles of nvr forward 
 *               products, forwardlmi[k] has space for deg+1 doubles;
 *   forwardlrg  is work space for the second lowest doubles of nvr forward 
 *               products, forwardlrg[k] has space for deg+1 doubles;
 *   forwardlpk  is work space for the lowest doubles of nvr forward 
 *               products, forwardlpk[k] has space for deg+1 doubles;
 *   backwardrtb is work space for the highest doubles of nvr-2 backward
 *               products, backwardrtb[k] can hold deg+1 doubles;
 *   backwardrix is work space for the second highest doubles of nvr-2
 *               backward products, backwardrix[k] can hold deg+1 doubles;
 *   backwardrmi is work space for the third highest doubles of nvr-2 
 *               backward products, backwardrmi[k] can hold deg+1 doubles;
 *   backwardrrg is work space for the fourth highest doubles of nvr-2
 *               backward products, backwardrrg[k] can hold deg+1 doubles;
 *   backwardrpk is work space for the fifth highest doubles of nvr-2
 *               backward products, backwardrpk[k] can hold deg+1 doubles;
 *   backwardltb is work space for the fifth lowest doubles of nvr-2
 *               backward products, backwardltb[k] can hold deg+1 doubles;
 *   backwardlix is work space for the fourth lowest doubles of nvr-2
 *               backward products, backwardlix[k] can hold deg+1 doubles;
 *   backwardlmi is work space for the third lowest doubles of nvr-2
 *               backward products, backwardlmi[k] can hold deg+1 doubles;
 *   backwardlrg is work space for the second lowest doubles of nvr-2
 *               backward products, backwardrg[k] can hold deg+1 doubles;
 *   backwardlpk is work space for the lowest doubles of nvr-2 backward
 *               products, backwardlpk[k] can hold deg+1 doubles;
 *   crossrtb    is work space for the highest doubles of nvr-2 cross
 *               products, crossrtb[k] has space for deg+1 doubles;
 *   crossrix    is work space for the second highest doubles of nvr-2
 *               cross products, crossrix[k] has space for deg+1 doubles;
 *   crossrmi    is work space for the third highest doubles of nvr-2
 *               cross products, crossrmi[k] has space for deg+1 doubles;
 *   crossrrg    is work space for the fourth highest doubles of nvr-2
 *               cross products, crossrrg[k] has space for deg+1 doubles;
 *   crossrpk    is work space for the fifth highest doubles of nvr-2
 *               cross products, crossrpk[k] has space for deg+1 doubles.
 *   crossltb    is work space for the fifth lowest doubles of nvr-2
 *               cross products, crossltb[k] has space for deg+1 doubles;
 *   crosslix    is work space for the fourth lowest doubles of nvr-2
 *               cross products, crosslix[k] has space for deg+1 doubles;
 *   crosslmi    is work space for the third lowest doubles of nvr-2
 *               cross products, crosslmi[k] has space for deg+1 doubles;
 *   crosslrg    is work space for the second lowest doubles of nvr-2
 *               cross products, crosslrg[k] has space for deg+1 doubles;
 *   crosslpk    is work space for the lowest doubles of nvr-2 cross
 *               products, crosslpk[k] has space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrtb  stores the highest doubles of the forward products,
 *   forwardrix  stores the second highest doubles of the forward products,
 *   forwardrmi  stores the third highest doubles of the forward products,
 *   forwardrrg  stores the fourth highest doubles of the forward products,
 *   forwardrpk  stores the fifth lowest doubles of the forward products,
 *   forwardltb  stores the fifth lowest doubles of the forward products,
 *   forwardlix  stores the fourth lowest doubles of the forward products,
 *   forwardlmi  stores the third lowest doubles of the forward products,
 *   forwardlrg  stores the second lowest doubles of the forward products,
 *   forwardlpk  stores the lowest doubles of the forward products,
 *               forward[nvr-1] contains the value of the product,
 *               forward[nvr-2] contains the derivative with respect
 *               to the last variable idx[nvr-1];
 *   backwardrtb stores the highest doubles of the backward products,
 *   backwardrix stores the second highest doubles of the backward products,
 *   backwardrmi stores the third highest doubles of the backward products,
 *   backwardrrg stores the fourth highest doubles of the backward products,
 *   backwardrpk stores the fifth highest doubles of the backward products,
 *   backwardltb stores the fifth lowest doubles of the backward products,
 *   backwardlix stores the fourth lowest doubles of the backward products,
 *   backwardlmi stores the third lowest doubles of the backward products,
 *   backwardlrg stores the second lowest doubles of the backward products,
 *   backwardlpk stores the lowest doubles of the backward products,
 *               backward[nvr-3] contains the derivative with respect
 *               to the first variable idx[0];
 *   crossrtb    stores the highest doubles of the cross products,
 *   crossrix    stores the second highest doubles of the cross products,
 *   crossrmi    stores the third highest doubles of the cross products,
 *   crossrrg    stores the fourth highest doubles of the cross products,
 *   crossrpk    stores the fifth highest doubles of the cross products,
 *   crossltb    stores the fifth lowest doubles of the cross products,
 *   crosslix    stores the fourth lowest doubles of the cross products,
 *   crosslmi    stores the third lowest doubles of the cross products,
 *   crosslrg    stores the second lowest doubles of the cross products,
 *   crosslpk    stores the lowest doubles of the cross products,
 *               cross[k] contains the derivatve with respect to
 *               variable idx[k+1]. */

void CPU_cmplx10_speel
 ( int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double **inputrertb, double **inputrerix, double **inputrermi,
   double **inputrerrg, double **inputrerpk,
   double **inputreltb, double **inputrelix, double **inputrelmi,
   double **inputrelrg, double **inputrelpk,
   double **inputimrtb, double **inputimrix, double **inputimrmi,
   double **inputimrrg, double **inputimrpk,
   double **inputimltb, double **inputimlix, double **inputimlmi,
   double **inputimlrg, double **inputimlpk,
   double **forwardrertb, double **forwardrerix, double **forwardrermi,
   double **forwardrerrg, double **forwardrerpk,
   double **forwardreltb, double **forwardrelix, double **forwardrelmi,
   double **forwardrelrg, double **forwardrelpk,
   double **forwardimrtb, double **forwardimrix, double **forwardimrmi,
   double **forwardimrrg, double **forwardimrpk,
   double **forwardimltb, double **forwardimlix, double **forwardimlmi,
   double **forwardimlrg, double **forwardimlpk,
   double **backwardrertb, double **backwardrerix, double **backwardrermi,
   double **backwardrerrg, double **backwardrerpk,
   double **backwardreltb, double **backwardrelix, double **backwardrelmi,
   double **backwardrelrg, double **backwardrelpk,
   double **backwardimrtb, double **backwardimrix, double **backwardimrmi,
   double **backwardimrrg, double **backwardimrpk,
   double **backwardimltb, double **backwardimlix, double **backwardimlmi,
   double **backwardimlrg, double **backwardimlpk,
   double **crossrertb, double **crossrerix, double **crossrermi,
   double **crossrerrg, double **crossrerpk,
   double **crossreltb, double **crossrelix, double **crossrelmi,
   double **crossrelrg, double **crossrelpk,
   double **crossimrtb, double **crossimrix, double **crossimrmi,
   double **crossimrrg, double **crossimrpk,
   double **crossimltb, double **crossimlix, double **crossimlmi,
   double **crossimlrg, double **crossimlpk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in penta double precision.
 *
 * REQUIRED : nvr > 2.
 *
 * ON ENTRY :
 *   nvr           number of variables in the product;
 *   deg           truncation degree of the series;
 *   idx           as many indices as the value of nvr,
 *                 idx[k] defines the place of the k-th variable,
 *                 with input values in input[idx[k]];
 *   cffrertb      highest doubles of the real parts of the coefficients
 *                 of the series of the product;
 *   cffrerix      second highest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrermi      third highest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrerrg      fourth highest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrerpk      fifth highest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffreltb      fifth lowest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrelix      fourth lowest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrelmi      third lowest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrelrg      second lowest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffrelpk      lowest doubles of the real parts of the
 *                 coefficients of the series of the product;
 *   cffimrtb      highest doubles of the imaginary pars of the
 *                 coefficients of the series of the product;
 *   cffimrix      second highest doubles of the imaginary pars of the
 *                 coefficients of the series of the product;
 *   cffimrmi      third highest doubles of the imaginary pars of the
 *                 coefficients of the series of the product;
 *   cffimrrg      fourth highest doubles of the imaginary pars of the
 *                 coefficients of the series of the product;
 *   cffimrpk      fifth highest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   cffimltb      fifth lowest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   cffimlix      fourth lowest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   cffimlmi      third lowest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   cffimlrg      second lowest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   cffimlpk      lowest doubles of the imaginary parts of the
 *                 coefficients of the series of the product;
 *   inputrertb    holds the highest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrerix    holds the second highest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrermi    holds the third highest of the real parts of
 *                 the coefficients of the input series;
 *   inputrerrg    holds the fourth highest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrerpk    holds the fifth highest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputreltb    holds the fifth lowest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrelix    holds the fourth lowest doubles of the real parts of 
 *                 the coefficients of the input series;
 *   inputrelmi    holds the third lowest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrelrg    holds the second lowest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputrelpk    holds the lowest doubles of the real parts of
 *                 the coefficients of the input series;
 *   inputimrtb    holds the highest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimrix    holds the second highest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimrmi    holds the third highest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimrrg    holds the fourth highest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimrpk    holds the fifth highest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimltb    holds the fifth lowest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimlix    holds the fourth lowest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimlmi    holds the third lowest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimlpk    holds the second lowest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   inputimlpk    holds the lowest doubles of the imaginary parts of
 *                 the coefficients of the input series;
 *   forwardrertb  is space for the highest doubles of nvr forward
 *                 products, forwardrertb[k] has space for deg+1 doubles;
 *   forwardrerix  is space for the second highest doubles of nvr forward
 *                 products, forwardrerix[k] has space for deg+1 doubles;
 *   forwardrermi  is space for the third highest doubles of nvr forward
 *                 products, forwardrermi[k] has space for deg+1 doubles;
 *   forwardrerrg  is space for the fourth highest doubles of nvr forward
 *                 products, forwardrerrg[k] has space for deg+1 doubles;
 *   forwardrerpk  is space for the fifth highest doubles of nvr forward
 *                 products, forwardrerpk[k] has space for deg+1 doubles;
 *   forwardreltb  is space for the fifth lowest doubles of nvr forward
 *                 products, forwardreltb[k] has space for deg+1 doubles;
 *   forwardrelix  is space for the fourth lowest doubles of nvr forward
 *                 products, forwardrelix[k] has space for deg+1 doubles;
 *   forwardrelmi  is space for the third lowest doubles of nvr forward
 *                 products, forwardrelmi[k] has space for deg+1 doubles;
 *   forwardrelrg  is space for the second lowest doubles of nvr forward
 *                 products, forwardrelrg[k] has space for deg+1 doubles;
 *   forwardrelpk  is space for the lowest doubles of nvr forward
 *                 products, forwardrelpk[k] has space for deg+1 doubles;
 *   forwardimltb  is space for the highest doubles of nvr forward
 *                 products, forwardimltb[k] has space for deg+1 doubles;
 *   forwardimlix  is space for the second highest doubles of nvr forward
 *                 products, forwardimlix[k] has space for deg+1 doubles;
 *   forwardimlmi  is space for the third highest doubles of nvr forward
 *                 products, forwardimlmi[k] has space for deg+1 doubles;
 *   forwardimlrg  is space for the fourth highest doubles of nvr forward
 *                 products, forwardimlrg[k] has space for deg+1 doubles;
 *   forwardimlpk  is space for the fifth highest doubles of nvr forward
 *                 products, forwardimlpk[k] has space for deg+1 doubles;
 *   forwardimrtb  is space for the fifth lowest doubles of nvr forward
 *                 products, forwardimrtb[k] has space for deg+1 doubles;
 *   forwardimrix  is space for the fourth lowest doubles of nvr forward
 *                 products, forwardimrix[k] has space for deg+1 doubles;
 *   forwardimrmi  is space for the third lowest doubles of nvr forward
 *                 products, forwardimrmi[k] has space for deg+1 doubles;
 *   forwardimrrg  is space for the second lowest doubles of nvr forward
 *                 products, forwardimrrg[k] has space for deg+1 doubles;
 *   forwardimrpk  is space for the lowest doubles of nvr forward
 *                 products, forwardimrpk[k] has space for deg+1 doubles;
 *   backwardrertb is space for the highest doubles of nvr-2 backward
 *                 products, backwardrertb[k] has space for deg+1 doubles;
 *   backwardrerix is space for the second highest doubles of nvr-2 backward
 *                 products, backwardrerix[k] has space for deg+1 doubles;
 *   backwardrermi is space for the third highest doubles of nvr-2 backward
 *                 products, backwardrermi[k] has space for deg+1 doubles;
 *   backwardrerrg is space for the fourth highest doubles of nvr-2 backward
 *                 products, backwardrerrg[k] has space for deg+1 doubles;
 *   backwardrerpk is space for the fifth highest doubles of nvr-2 backward
 *                 products, backwardrerpk[k] has space for deg+1 doubles;
 *   backwardreltb is space for the fifth lowest doubles of nvr-2 backward
 *                 products, backwardreltb[k] has space for deg+1 doubles;
 *   backwardrelix is space for the fourth lowest doubles of nvr-2 backward
 *                 products, backwardrelix[k] has space for deg+1 doubles;
 *   backwardrelmi is space for the third lowest doubles of nvr-2 backward
 *                 products, backwardrelmi[k] has space for deg+1 doubles;
 *   backwardrelrg is space for the second lowest doubles of nvr-2 backward
 *                 products, backwardrelrg[k] has space for deg+1 doubles;
 *   backwardrelpk is space for the lowest doubles of nvr-2 backward
 *                 products, backwardrelpk[k] has space for deg+1 doubles;
 *   backwardimrtb is space for the highest doubles of nvr-2 backward
 *                 products, backwardimrtb[k] has space for deg+1 doubles;
 *   backwardimrix is space for the second highest doubles of nvr-2 backward
 *                 products, backwardimrix[k] has space for deg+1 doubles;
 *   backwardimrmi is space for the third highest doubles of nvr-2 backward
 *                 products, backwardimrmi[k] has space for deg+1 doubles;
 *   backwardimrrg is space for the fourth highest doubles of nvr-2 backward
 *                 products, backwardimrrg[k] has space for deg+1 doubles;
 *   backwardimrpk is space for the fifth highest doubles of nvr-2 backward
 *                 products, backwardimrpk[k] has space for deg+1 doubles;
 *   backwardimltb is space for the fifth lowest doubles of nvr-2 backward
 *                 products, backwardimltb[k] has space for deg+1 doubles;
 *   backwardimlix is space for the fourth lowest doubles of nvr-2 backward
 *                 products, backwardimlix[k] has space for deg+1 doubles;
 *   backwardimlmi is space for the third lowest doubles of nvr-2 backward
 *                 products, backwardimlmi[k] has space for deg+1 doubles;
 *   backwardimlrg is space for the second lowest doubles of nvr-2 backward
 *                 products, backwardimlrg[k] has space for deg+1 doubles;
 *   backwardimlpk is space for the lowest doubles of nvr-2 backward
 *                 products, backwardimlpk[k] has space for deg+1 doubles;
 *   crossrertb    is space for the highest doubles of nvr-2 cross
 *                 products, crossretb[k] has space for deg+1 doubles;
 *   crossrerix    is space for the second highest doubles of nvr-2 cross
 *                 products, crossreix[k] has space for deg+1 doubles;
 *   crossrermi    is space for the third highest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossrerrg    is space for the fourth highest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossrerpk    is space for the fifth highest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossreltb    is space for the fifth lowest doubles of nvr-2 cross
 *                 products, crossretb[k] has space for deg+1 doubles;
 *   crossrelix    is space for the fourth lowest doubles of nvr-2 cross
 *                 products, crossreix[k] has space for deg+1 doubles;
 *   crossrelmi    is space for the third lowest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossrelrg    is space for the second lowest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossrelpk    is space for the lowest doubles of nvr-2 cross
 *                 products, crossrehi[k] has space for deg+1 doubles;
 *   crossimrtb    is space for the highest doubles of nvr-2 cross
 *                 products, crossimrtb[k] has space for deg+1 doubles;
 *   crossimrix    is space for the second highest doubles of nvr-2 cross
 *                 products, crossimrix[k] has space for deg+1 doubles;
 *   crossimrmi    is space for the third highest doubles of nvr-2 cross
 *                 products, crossimrmi[k] has space for deg+1 doubles;
 *   crossimrrg    is space for the fourth highest doubles of nvr-2 cross
 *                 products, crossimrrg[k] has space for deg+1 doubles.
 *   crossimrpk    is space for the fifth highest doubles of nvr-2 cross
 *                 products, crossimrpk[k] has space for deg+1 doubles;
 *   crossimltb    is space for the fifth lowest doubles of nvr-2 cross
 *                 products, crossimltb[k] has space for deg+1 doubles;
 *   crossimlix    is space for the fourth lowest doubles of nvr-2 cross
 *                 products, crossimlix[k] has space for deg+1 doubles;
 *   crossimlmi    is space for the third lowest doubles of nvr-2 cross
 *                 products, crossimlmi[k] has space for deg+1 doubles;
 *   crossimlrg    is space for the second lowest doubles of nvr-2 cross
 *                 products, crossimlrg[k] has space for deg+1 doubles.
 *   crossimlpk    is space for the lowest doubles of nvr-2 cross
 *                 products, crossimlpk[k] has space for deg+1 doubles.
 *
 * ON RETURN :
 *   forwardrertb  stores the highest doubles of the real parts
 *                 of the forward products,
 *   forwardrerix  stores the second highest doubles of the real parts
 *                 of the forward products,
 *   forwardrermi  stores the third highest doubles of the real parts
 *                 of the forward products,
 *   forwardrerrg  stores the fourth highest doubles of the real parts
 *                 of the forward products,
 *   forwardrerpk  stores the fifth highest doubles of the real parts
 *                 of the forward products,
 *   forwardreltb  stores the fifth lowest doubles of the real parts
 *                 of the forward products,
 *   forwardrelix  stores the fourth lowest doubles of the real parts
 *                 of the forward products,
 *   forwardrelmi  stores the third lowest doubles of the real parts
 *                 of the forward products,
 *   forwardrelrg  stores the second lowest doubles of the real parts
 *                 of the forward products,
 *   forwardrelpk  stores the lowest doubles of the real parts
 *                 of the forward products,
 *   forwardimrtb  stores the highest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimrix  stores the second highest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimrmi  stores the third highest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimrrg  stores the fourth highest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimrpk  stores the fifth highest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimltb  stores the fifth lowest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimlix  stores the fourth lowest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimlmi  stores the third lowest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimlrg  stores the second lowest doubles of the imaginary parts
 *                 of the forward products,
 *   forwardimlpk  stores the lowest doubles of the imaginary parts
 *                 of the forward products,
 *                 forward[nvr-1] contains the value of the product,
 *                 forward[nvr-2] contains the derivative with respect
 *                 to the last variable idx[nvr-1];
 *   backwardrertb stores the highest doubles of the real parts
 *                 of the backward products,
 *   backwardrerix stores the second highest doubles of the real parts
 *                 of the backward products,
 *   backwardrermi stores the third highest doubles of the real parts
 *                 of the backward products,
 *   backwardrerrg stores the fourth highest doubles of the real parts
 *                 of the backward products,
 *   backwardrerpk stores the fifth highest doubles of the real parts
 *                 of the backward products,
 *   backwardreltb stores the fifth lowest doubles of the real parts
 *                 of the backward products,
 *   backwardrelix stores the fourth lowest doubles of the real parts
 *                 of the backward products,
 *   backwardrelmi stores the third lowest doubles of the real parts
 *                 of the backward products,
 *   backwardrelrg stores the second lowest doubles of the real parts
 *                 of the backward products,
 *   backwardrelpk stores the lowest doubles of the real parts
 *                 of the backward products,
 *   backwardimrtb stores the highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimrix stores the second highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimrmi stores the third highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimrrg stores the fourth highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimrpk stores the fifth highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimltb stores the fifth lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimlix stores the fourth lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimlmi stores the third lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimlrg stores the second lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardimlpk stores the lowest doubles of the imaginary parts
 *                 of the backward products,
 *                 backward[nvr-3] contains the derivative with respect
 *                 to the first variable idx[0];
 *   crossrertb    stores the highest doubles of the real parts
 *                 of the cross products,
 *   crossrerix    stores the second highest doubles of the real parts
 *                 of the cross products,
 *   crossrermi    stores the third highest doubles of the real parts
 *                 of the cross products,
 *   crossrerrg    stores the fourth highest doubles of the real parts
 *                 of the cross products,
 *   crossrerpk    stores the fifth highest doubles of the real parts
 *                 of the cross products,
 *   crossreltb    stores the fifth lowest doubles of the real parts
 *                 of the cross products,
 *   crossrelix    stores the fourth lowest doubles of the real parts
 *                 of the cross products,
 *   crossrelmi    stores the third lowest doubles of the real parts
 *                 of the cross products,
 *   crossrelrg    stores the second lowest doubles of the real parts
 *                 of the cross products,
 *   crossrelpk    stores the lowest doubles of the real parts
 *                 of the cross products,
 *   crossimrtb    stores the highest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimrix    stores the second highest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimrmi    stores the third highest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimrrg    stores the fourth highest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimrpk    stores the fifth highest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimltb    stores the fifth lowest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimlix    stores the fourth lowest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimlmi    stores the third lowest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimlrg    stores the second lowest doubles of the imaginary parts
 *                 of the cross products,
 *   crossimlpk    stores the lowest doubles of the imaginary parts
 *                 of the cross products,
 *                 cross[k] contains the derivatve with respect to
 *                 the variable idx[k+1]. */

void CPU_dbl10_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   dim       total number of variables;
 *   deg       truncation degree of the series;
 *   idx       as many indices as the value of nvr,
 *             idx[k] defines the place of the k-th variable,
 *             with input values in input[idx[k]];
 *   cffrtb    deg+1 higest doubles in the coefficient series;
 *   cffrix    deg+1 second higest doubles in the coefficient series;
 *   cffrmi    deg+1 third highest doubles in the coefficient series;
 *   cffrrg    deg+1 fourth highest doubles in the coefficient series;
 *   cffrpk    deg+1 fifth highest doubles in the coefficient series;
 *   cffltb    deg+1 fifth lowest doubles in the coefficient series;
 *   cfflix    deg+1 fourth lowest doubles in the coefficient series;
 *   cfflmi    deg+1 third lowest doubles in the coefficient series;
 *   cfflrg    deg+1 second lowest doubles in the coefficient series;
 *   cfflpk    deg+1 lowest doubles in the coefficient series;
 *   inputrtb  holds the highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrix  holds the second highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrmi  holds the third highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrrg  holds the fourth highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrpk  holds the fifth highest doubles of the input series
 *             for all variables in the monomial;
 *   inputltb  holds the fifth lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlix  holds the fourth lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlmi  holds the third lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlrg  holds the second lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlpk  holds the lowest doubles of the input series
 *             for all variables in the monomial;
 *   outputrtb has space allocated for dim+1 series of degree deg;
 *   outputrix has space allocated for dim+1 series of degree deg;
 *   outputrmi has space allocated for dim+1 series of degree deg;
 *   outputrrg has space allocated for dim+1 series of degree deg;
 *   outputrpk has space allocated for dim+1 series of degree deg;
 *   outputltb has space allocated for dim+1 series of degree deg;
 *   outputlix has space allocated for dim+1 series of degree deg;
 *   outputlmi has space allocated for dim+1 series of degree deg;
 *   outputlrg has space allocated for dim+1 series of degree deg;
 *   outputlpk has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputrtb stores highest doubles of the derivatives and the value,
 *   outputrix stores second highest doubles of the derivatives and the value,
 *   outputrmi stores third highest doubles of the derivatives and the value,
 *   outputrrg stores fourth highest doubles of the derivatives and the value,
 *   outputrpk stores fifth highest doubles of the derivatives and the value,
 *   outputltb stores fifth lowest doubles of the derivatives and the value,
 *   outputlix stores fourth lowest doubles of the derivatives and the value,
 *   outputlmi stores third lowest doubles of the derivatives and the value,
 *   outputlrg stores second lowest doubles of the derivatives and the value,
 *   outputlpk stores lowest doubles of the derivatives and the value,
 *             output[idx[k]], for k from 0 to nvr, contains the
 *             deriviative with respect to the variable idx[k];
 *             output[dim] contains the value of the product. */

void CPU_cmplx10_evaldiff
 ( int dim, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double **inputrertb, double **inputrerix, double **inputrermi,
   double **inputrerrg, double **inputrerpk,
   double **inputreltb, double **inputrelix, double **inputrelmi,
   double **inputrelrg, double **inputrelpk,
   double **inputimrtb, double **inputimrix, double **inputimrmi,
   double **inputimrrg, double **inputimrpk,
   double **inputimltb, double **inputimlix, double **inputimlmi,
   double **inputimlrg, double **inputimlpk,
   double **outputrertb, double **outputrerix, double **outputrermi,
   double **outputrerrg, double **outputrerpk,
   double **outputreltb, double **outputrelix, double **outputrelmi,
   double **outputrelrg, double **outputrelpk,
   double **outputimrtb, double **outputimrix, double **outputimrmi,
   double **outputimrrg, double **outputimrpk,
   double **outputimltb, double **outputimlix, double **outputimlmi,
   double **outputimlrg, double **outputimlpk );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a monomial.
 *
 * ON ENTRY :
 *   dim         total number of variables;
 *   deg         truncation degree of the series;
 *   idx         as many indices as the value of nvr,
 *               idx[k] defines the place of the k-th variable,
 *               with input values in input[idx[k]];
 *   cffrertb    highest doubles of the real parts of the coefficients
 *               of the series of the product;
 *   cffrerix    second highest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrermi    third highest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrerrg    fourth highest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrerpk    fifth highest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffreltb    fifth lowest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrelix    fourth lowest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrelmi    third lowest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrelrg    second lowest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffrelpk    lowest doubles of the real parts of the
 *               coefficients of the series of the product;
 *   cffimrtb    highest doubles of the imaginary pars of the
 *               coefficients of the series of the product;
 *   cffimrix    second highest doubles of the imaginary pars of the
 *               coefficients of the series of the product;
 *   cffimrmi    third highest doubles of the imaginary pars of the
 *               coefficients of the series of the product;
 *   cffimrrg    fourth highest doubles of the imaginary pars of the
 *               coefficients of the series of the product;
 *   cffimrpk    fifth highest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   cffimltb    fifth lowest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   cffimlix    fourth lowest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   cffimlmi    third lowest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   cffimlrg    second lowest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   cffimlpk    lowest doubles of the imaginary parts of the
 *               coefficients of the series of the product;
 *   inputrertb  holds the highest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrerix  holds the second highest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrermi  holds the third highest of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrerrg  holds the fourth highest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrerpk  holds the fifth highest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputreltb  holds the fifth lowest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrelix  holds the fourth lowest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrelmi  holds the third lowest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrelrg  holds the second lowest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputrelpk  holds the lowest doubles of the real parts of the 
 *               coefficients of the series for all variables in the monomial;
 *   inputimrtb  holds the highest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimrix  holds the second highest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimrmi  holds the third highest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimrrg  holds the fourth highest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimrpk  holds the fifth highest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimltb  holds the fifth lowest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimlix  holds the fourth lowest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimlmi  holds the third lowest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimlpk  holds the second lowest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   inputimlpk  holds the lowest doubles of the imaginary parts of the
 *               coefficients of the series for all variables in the monomial;
 *   outputrertb has space allocated for dim+1 series of degree deg;
 *   outputrerix has space allocated for dim+1 series of degree deg;
 *   outputrermi has space allocated for dim+1 series of degree deg;
 *   outputrerrg has space allocated for dim+1 series of degree deg;
 *   outputrerpk has space allocated for dim+1 series of degree deg;
 *   outputreltb has space allocated for dim+1 series of degree deg;
 *   outputrelix has space allocated for dim+1 series of degree deg;
 *   outputrelmi has space allocated for dim+1 series of degree deg;
 *   outputrelrg has space allocated for dim+1 series of degree deg;
 *   outputrelpk has space allocated for dim+1 series of degree deg;
 *   outputimrtb has space allocated for dim+1 series of degree deg;
 *   outputimrix has space allocated for dim+1 series of degree deg;
 *   outputimrmi has space allocated for dim+1 series of degree deg;
 *   outputimrrg has space allocated for dim+1 series of degree deg;
 *   outputimrpk has space allocated for dim+1 series of degree deg;
 *   outputimltb has space allocated for dim+1 series of degree deg;
 *   outputimlix has space allocated for dim+1 series of degree deg;
 *   outputimlmi has space allocated for dim+1 series of degree deg;
 *   outputimlrg has space allocated for dim+1 series of degree deg;
 *   outputimlpk has space allocated for dim+1 series of degree deg.
 *
 * ON RETURN :
 *   outputrertb stores the highest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrerix stores the second highest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrermi stores the third highest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrerrg stores the fourth highest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrerpk stores the fifth highest doubles of the real parts
 *               of the derivatives and the value,
 *   outputreltb stores the fifth lowest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrelix stores the fourth lowest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrelmi stores the third lowest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrelrg stores the second lowest doubles of the real parts
 *               of the derivatives and the value,
 *   outputrelpk stores the lowest doubles of the real parts
 *               of the derivatives and the value,
 *   outputimrtb stores the highest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimrix stores the second highest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimrmi stores the third highest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimrrg stores the fourth highest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimrpk stores the fifth highest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimltb stores the fifth lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlix stores the fourth lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlmi stores the third lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlrg stores the second lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlpk stores the lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *               output[idx[k]], for k from 0 to nvr, contains the
 *               deriviative with respect to the variable idx[k];
 *               output[dim] contains the value of the product. */

#endif
