// The file dbl10_monomials_kernels.h specifies functions to evaluate and
// differentiate a monomial at power series truncated to the same degree,
// in deca double precision.

#ifndef __dbl10_monomials_kernels_h__
#define __dbl10_monomials_kernels_h__

void GPU_dbl10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double *inputrtb, double *inputrix, double *inputrmi, double *inputrrg,
   double *inputrpk, double *inputltb, double *inputlix, double *inputlmi,
   double *inputlrg, double *inputlpk,
   double *forwardrtb, double *forwardrix, double *forwardrmi,
   double *forwardrrg, double *forwardrpk,
   double *forwardltb, double *forwardlix, double *forwardlmi,
   double *forwardlrg, double *forwardlpk,
   double *backwardrtb, double *backwardrix, double *backwardrmi,
   double *backwardrrg, double *backwardrpk,
   double *backwardltb, double *backwardlix, double *backwardlmi,
   double *backwardlrg, double *backwardlpk,
   double *crossrtb, double *crossrix, double *crossrmi, double *crossrrg,
   double *crossrpk, double *crossltb, double *crosslix, double *crosslmi,
   double *crosslrg, double *crosslpk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for real coefficients in deca double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS          number of threads in one block, must be deg+1;
 *   nvr         number of variables in the product;
 *   deg         truncation degree of the series;
 *   idx         as many indices as the value of nvr,
 *               idx[k] defines the place of the k-th variable,
 *               with input values in input[idx[k]];
 *   cffrtb      deg+1 highest coefficient doubles;
 *   cffrix      deg+1 second highest coefficient doubles;
 *   cffrmi      deg+1 third highest coefficient doubles;
 *   cffrrg      deg+1 fourth highest coefficient doubles;
 *   cffrpk      deg+1 fifth highest coefficient doubles;
 *   cffltb      deg+1 fifth lowest coefficient doubles;
 *   cfflix      deg+1 fourth lowest coefficient doubles;
 *   cfflmi      deg+1 third lowest coefficient doubles;
 *   cfflrg      deg+1 second lowest coefficient doubles;
 *   cfflpk      deg+1 lowest coefficient doubles.
 *   inputrtb    stores the highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrix    stores the second highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrmi    stores the third highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrrg    stores the fourth highest doubles of the input series
 *               for all variables in the monomial;
 *   inputrpk    stores the fifth highest doubles of the input series
 *               for all variables in the monomial;
 *   inputltb    stores the fifth lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlix    stores the fourth lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlmi    stores the third lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlrg    stores the second lowest doubles of the input series
 *               for all variables in the monomial;
 *   inputlpk    stores the lowest doubles of the input series
 *               for all variables in the monomial;
 *   forwardrtb  is work space for the highest doubles of nvr 
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardrix  is work space for the second highest doubles of nvr 
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardrmi  is work space for the third highest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardrrg  is work space for the fourth highest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardrpk  is work space for the fifth highest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardltb  is work space for the fifth lowest doubles of nvr 
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardlix  is work space for the fourth lowest doubles of nvr 
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardlmi  is work space for the third lowest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardlrg  is work space for the second lowest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   forwardlpk  is work space for the lowest doubles of nvr
 *               forward products, has space for nvr*(deg+1) doubles;
 *   backwardrtb is work space for the highest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardrix is work space for the second highest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardrmi is work space for the third highest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardrrg is work space for the fourth highest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardrpk is work space for the fifth highest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardltb is work space for the fifth lowest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardlix is work space for the fourth lowest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardlmi is work space for the third lowest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardlrg is work space for the second lowest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   backwardlpk is work space for the lowest doubles of nvr-2
 *               backward products, has space for (nvr-2)*(deg+1) doubles;
 *   crossrtb    is work space for the highest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossrix    is work space for the second highest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossrmi    is work space for the third highest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crossrrg    is work space for the fourth highest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles.
 *   crossrpk    is work space for the fifth highest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles.
 *   crossltb    is work space for the fifth lowest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crosslix    is work space for the fourth lowest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crosslmi    is work space for the third lowest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles;
 *   crosslrg    is work space for the second lowest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles.
 *   crosslpk    is work space for the lowest doubles of nvr-2
 *               cross products, has space for (nvr-2)*(deg+1) doubles.
 *
 * ON RETURN :
 *   forwardrtb  stores the highest doubles of the forward products,
 *   forwardrix  stores the second highest doubles of the forward products,
 *   forwardrmi  stores the third highest doubles of the forward products,
 *   forwardrrg  stores the fourth highest doubles of the forward products,
 *   forwardrpk  stores the fifth highest doubles of the forward products,
 *   forwardltb  stores the fifth lowest doubles of the forward products,
 *   forwardlix  stores the fourth lowest doubles of the forward products,
 *   forwardlmi  stores the third lowest doubles of the forward products,
 *   forwardlrg  stores the second lowest doubles of the forward products,
 *   forwardlpk  stores the lowest doubles of the forward products,
 *               forward[nvr-1] contains the value of the product,
 *               forward[nvr-2] contains the derivative with respect
 *               to the last variable idx[nvr-1] if nvr > 2;
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
 *               to the first variable idx[0] if nvr > 2;
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
 *               cross[k] contains the derivative with respect to
 *               variable idx[k+1]. */

void GPU_cmplx10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double *inputrertb, double *inputrerix, double *inputrermi,
   double *inputrerrg, double *inputrerpk,
   double *inputreltb, double *inputrelix, double *inputrelmi,
   double *inputrelrg, double *inputrelpk,
   double *inputimrtb, double *inputimrix, double *inputimrmi,
   double *inputimrrg, double *inputimrpk,
   double *inputimltb, double *inputimlix, double *inputimlmi,
   double *inputimlrg, double *inputimlpk,
   double *forwardrertb, double *forwardrerix, double *forwardrermi,
   double *forwardrerrg, double *forwardrerpk,
   double *forwardreltb, double *forwardrelix, double *forwardrelmi,
   double *forwardrelrg, double *forwardrelpk,
   double *forwardimrtb, double *forwardimrix, double *forwardimrmi,
   double *forwardimrrg, double *forwardimrpk,
   double *forwardimltb, double *forwardimlix, double *forwardimlmi,
   double *forwardimlrg, double *forwardimlpk,
   double *backwardrertb, double *backwardrerix, double *backwardrermi,
   double *backwardrerrg, double *backwardrerpk,
   double *backwardreltb, double *backwardrelix, double *backwardrelmi,
   double *backwardrelrg, double *backwardrelpk,
   double *backwardimrtb, double *backwardimrix, double *backwardimrmi,
   double *backwardimrrg, double *backwardimrpk,
   double *backwardimltb, double *backwardimlix, double *backwardimlmi,
   double *backwardimlrg, double *backwardimlpk,
   double *crossrertb, double *crossrerix, double *crossrermi,
   double *crossrerrg, double *crossrerpk,
   double *crossreltb, double *crossrelix, double *crossrelmi,
   double *crossrelrg, double *crossrelpk,
   double *crossimrtb, double *crossimrix, double *crossimrmi,
   double *crossimrrg, double *crossimrpk,
   double *crossimltb, double *crossimlix, double *crossimlmi,
   double *crossimlrg, double *crossimlpk );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a product of variables at power series truncated to the same degree,
 *   for complex coefficients in deca double precision.
 *
 * REQUIRED : nvr > 2 and BS = deg+1.
 *   The cff and input are allocated on the device
 *   and all coefficients and input series are copied from host to device.
 *   The forward, backward, and cross are allocated on the device.
 *
 * ON ENTRY :
 *   BS            number of threads in one block, must be deg+1;
 *   nvr           number of variables in the product;
 *   deg           truncation degree of the series;
 *   idx           as many indices as the value of nvr,
 *                 idx[k] defines the place of the k-th variable,
 *                 with input values in input[idx[k]];
 *   cffrertb      highest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrerix      second highest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrermi      middle doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrerrg      second lowest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrerpk      lowest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffreltb      highest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrelix      second highest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrelmi      middle doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrelrg      second lowest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffrelpk      lowest doubles of the real parts
 *                 of the series coefficient of the product;
 *   cffimrtb      highest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimrix      second highest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimrmi      third highest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimrrg      fourth highest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimrpk      fifth highest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimltb      fifth lowest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimlix      fourth lowest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimlmi      third lowest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimlrg      second lowest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   cffimlpk      lowest doubles of the imaginary parts
 *                 of the series coefficient of the product;
 *   inputrertb    stores the highest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrerix    stores the second highest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrermi    stores the third highest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrerrg    stores the fourth highest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrerpk    stores the fifth highest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputreltb    stores the fifth lowest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrelix    stores the fourth lowest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrelmi    stores the third lowest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrelrg    stores the second lowest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputrelpk    stores the lowest doubles of the real parts
 *                 of the coefficients of the input series;
 *   inputimtb     stores the highest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimrix    stores the second highest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimrmi    stores the third highest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimrrg    stores the fourth highest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimrpk    stores the fifth highest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimltb    stores the fifth lowest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimlix    stores the fourth lowest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimlmi    stores the third lowest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimlrg    stores the second lowest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   inputimlpk    stores the lowest doubles of the imaginary parts
 *                 of the coefficients of the input series;
 *   forwardrertb  is work space for the highest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrertb has space for nvr*(deg+1) doubles;
 *   forwardrerix  is work space for the second highest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrerix has space for nvr*(deg+1) doubles;
 *   forwardrermi  is work space for the third highest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrermi has space for nvr*(deg+1) doubles;
 *   forwardrerrg  is work space for the fourth highest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrerrg has space for nvr*(deg+1) doubles;
 *   forwardrerpk  is work space for the fifth highest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrerpk has space for nvr*(deg+1) doubles;
 *   forwardreltb  is work space for the fifth lowest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardreltb has space for nvr*(deg+1) doubles;
 *   forwardrelix  is work space for the fourth lowest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrelix has space for nvr*(deg+1) doubles;
 *   forwardrelmi  is work space for the third lowest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrelmi has space for nvr*(deg+1) doubles;
 *   forwardrelrg  is work space for the second lowest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrelrg has space for nvr*(deg+1) doubles;
 *   forwardrelpk  is work space for the lowest doubles
 *                 for all real parts of the nvr forward products,
 *                 forwardrelpk has space for nvr*(deg+1) doubles;
 *   forwardimrtb  is work space for the highest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimrtb has space for nvr*(deg+1) doubles;
 *   forwardimrix  is work space for the second highest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimrix has space for nvr*(deg+1) doubles;
 *   forwardimrmi  is work space for the third highest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimrmi has space for nvr*(deg+1) doubles;
 *   forwardimrrg  is work space for the fourth highest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimrrg is space for nvr*(deg+1) doubles;
 *   forwardimrpk  is work space for the fifth highest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimrpk is space for nvr*(deg+1) doubles;
 *   forwardimltb  is work space for the fifth lowest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimltb has space for nvr*(deg+1) doubles;
 *   forwardimlix  is work space for the fourth lowest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimlix has space for nvr*(deg+1) doubles;
 *   forwardimlmi  is work space for the third lowest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimlmi has space for nvr*(deg+1) doubles;
 *   forwardimlrg  is work space for the second lowest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimlrg is space for nvr*(deg+1) doubles;
 *   forwardimlpk  is work space for the lowest doubles
 *                 for all imaginary parts of the nvr forward products,
 *                 forwardimlpk is space for nvr*(deg+1) doubles;
 *   backwardrertb is work space for all highest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrertb has space for (nvr-2)*(deg+1) doubles;
 *   backwardrerix is work space for all second highest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrerix has space for (nvr-2)*(deg+1) doubles;
 *   backwardrermi is work space for all third highest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrermi has space for (nvr-2)*(deg+1) doubles;
 *   backwardrerrg is work space for all fourth highest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrerrg has space for (nvr-2)*(deg+1) doubles;
 *   backwardrerpk is work space for all fifth highest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrerpk has space for (nvr-2)*(deg+1) doubles;
 *   backwardreltb is work space for all fifth lowest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardreltb has space for (nvr-2)*(deg+1) doubles;
 *   backwardrelix is work space for all fourth lowest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrelix has space for (nvr-2)*(deg+1) doubles;
 *   backwardrelmi is work space for all third lowest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrelmi has space for (nvr-2)*(deg+1) doubles;
 *   backwardrelrg is work space for all second lowest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrelrg has space for (nvr-2)*(deg+1) doubles;
 *   backwardrelpk is work space for all lowest doubles
 *                 for all real parts of the nvr-2 backward products,
 *                 backwardrelpk has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrtb is work space for the highest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimrtb has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrix is work space for the second highest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimrix has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrmi is work space for the third highest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimrmi has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrrg is work space for the fourth highest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimrrg has space for (nvr-2)*(deg+1) doubles;
 *   backwardimrpk is work space for the fifth highest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimrpk has space for (nvr-2)*(deg+1) doubles;
 *   backwardimltb is work space for the fifth lowest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimltb has space for (nvr-2)*(deg+1) doubles;
 *   backwardimlix is work space for the fourth lowest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimlix has space for (nvr-2)*(deg+1) doubles;
 *   backwardimlmi is work space for the third lowest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimlmi has space for (nvr-2)*(deg+1) doubles;
 *   backwardimlrg is work space for the second lowest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimlrg has space for (nvr-2)*(deg+1) doubles;
 *   backwardimlpk is work space for the lowest doubles
 *                 for all imaginary parts of the nvr-2 backward products,
 *                 backwardimlpk has space for (nvr-2)*(deg+1) doubles;
 *   crossrertb    is work space for the highest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrertb has space for (nvr-2)*(deg+1) doubles;
 *   crossrerix    is work space for the second highest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrerix has space for (nvr-2)*(deg+1) doubles;
 *   crossrermi    is work space for the third highest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrermi has space for (nvr-2)*(deg+1) doubles;
 *   crossrerrg    is work space for the fourth highest doubles
 *                 for the real parts of nvr-2 cross products,
 *                 crossrerrg has space for (nvr-2)*(deg+1) doubles;
 *   crossrerpk    is work space for the fifth highest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrerpk has space for (nvr-2)*(deg+1) doubles;
 *   crossreltb    is work space for the fifth lowest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossreltb has space for (nvr-2)*(deg+1) doubles;
 *   crossrelix    is work space for the fourth lowest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrelix has space for (nvr-2)*(deg+1) doubles;
 *   crossrelmi    is work space for the third lowest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrelmi has space for (nvr-2)*(deg+1) doubles;
 *   crossrelrg    is work space for the second lowest doubles
 *                 for the real parts of nvr-2 cross products,
 *                 crossrelrg has space for (nvr-2)*(deg+1) doubles;
 *   crossrelpk    is work space for the lowest doubles
 *                 for the real parts of the nvr-2 cross products,
 *                 crossrelpk has space for (nvr-2)*(deg+1) doubles;
 *   crossimrtb    is work space for the highest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimrtb has space for (nvr-2)*(deg+1) doubles;
 *   crossimrix    is work space for the second highest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimrix has space for (nvr-2)*(deg+1) doubles;
 *   crossimrmi    is work space for the third highest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimrmi has space for (nvr-2)*(deg+1) doubles;
 *   crossimrrg    is work space for the fourth highest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimrrg has space for (nvr-2)*(deg+1) doubles;
 *   crossimrpk    is work space for the fifth highest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimrpk has space for (nvr-2)*(deg+1) doubles.
 *   crossimltb    is work space for the fifth lowest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimltb has space for (nvr-2)*(deg+1) doubles;
 *   crossimlix    is work space for the fourth lowest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimlix has space for (nvr-2)*(deg+1) doubles;
 *   crossimlmi    is work space for the third lowest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimlmi has space for (nvr-2)*(deg+1) doubles;
 *   crossimlrg    is work space for the second lowest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimlrg has space for (nvr-2)*(deg+1) doubles;
 *   crossimlpk    is work space for the lowest doubles
 *                 for the imaginary parts of nvr-2 cross products,
 *                 crossimlpk has space for (nvr-2)*(deg+1) doubles.
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
 *   forwardreltb  stores the second lowest doubles of the real parts
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
 *                 of the forward products
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
 *   backwardrertb stores the highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrerix stores the second highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrermi stores the third highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrerrg stores the fourth highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrerpk stores the fifth highest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardreltb stores the fifth lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrelix stores the fourth lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrelmi stores the third lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrelrg stores the second lowest doubles of the imaginary parts
 *                 of the backward products,
 *   backwardrelpk stores the lowest doubles of the imaginary parts of 
 *                 the backward products,
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
 *                 cross[k] contains the derivative with respect to
 *                 variable idx[k+1]. */

void GPU_dbl10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
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
 *   BS        number of threads in one block, must be deg+1;
 *   dim       total number of variables;
 *   deg       truncation degree of the series;
 *   idx       as many indices as the value of nvr,
 *             idx[k] defines the place of the k-th variable,
 *             with input values in input[idx[k]];
 *   cffrtb    deg+1 highest coefficient doubles;
 *   cffrix    deg+1 second highest coefficient doubles;
 *   cffrmi    deg+1 third highest coefficient doubles;
 *   cffrrg    deg+1 fourth highest coefficient doubles;
 *   cffrpk    deg+1 fifth highest coefficient doubles;
 *   cffltb    deg+1 fifth lowest coefficient doubles;
 *   cfflix    deg+1 fourth lowest coefficient doubles;
 *   cfflmi    deg+1 third lowest coefficient doubles;
 *   cfflrg    deg+1 second lowest coefficient doubles;
 *   cfflpk    deg+1 lowest coefficient doubles.
 *   inputrtb  stores the highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrix  stores the second highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrmi  stores the third highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrrg  stores the fourth highest doubles of the input series
 *             for all variables in the monomial;
 *   inputrpk  stores the fifth highest doubles of the input series
 *             for all variables in the monomial;
 *   inputltb  stores the fifth lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlix  stores the fourth lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlmi  stores the third lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlrg  stores the second lowest doubles of the input series
 *             for all variables in the monomial;
 *   inputlpk  stores the lowest doubles of the input series
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
 *   outputrtb stores the highest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputrix stores the second highest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputrmi stores the third highest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputrrg stores the fourth highest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputrpk stores the fifth highest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputltb stores the fifth lowest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputlix stores the fourth lowest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputlmi stores the third lowest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputlrg stores the second lowest doubles of the derivatives and 
 *             the value of the monomial,
 *   outputlpk stores the lowest doubles of the derivatives and 
 *             the value of the monomial,
 *             output[idx[k]], for k from 0 to nvr, contains the
 *             deriviative with respect to the variable idx[k];
 *             output[dim] contains the value of the product. */

void GPU_cmplx10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi,
   double *cffrerrg, double *cffrerpk,
   double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi,
   double *cffimrrg, double *cffimrpk,
   double *cffimltb, double *cffimlix, double *cffimlmi,
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
 *   BS          number of threads in one block, must be deg+1;
 *   dim         total number of variables;
 *   deg         truncation degree of the series;
 *   idx         as many indices as the value of nvr,
 *               idx[k] defines the place of the k-th variable,
 *               with input values in input[idx[k]];
 *   cffrertb    highest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrerix    second highest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrermi    middle doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrerrg    second lowest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrerpk    lowest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffreltb    highest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrelix    second highest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrelmi    middle doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrelrg    second lowest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffrelpk    lowest doubles of the real parts
 *               of the series coefficient of the product;
 *   cffimrtb    highest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimrix    second highest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimrmi    third highest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimrrg    fourth highest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimrpk    fifth highest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimltb    fifth lowest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimlix    fourth lowest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimlmi    third lowest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimlrg    second lowest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   cffimlpk    lowest doubles of the imaginary parts
 *               of the series coefficient of the product;
 *   inputrertb  stores the highest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrerix  stores the second highest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrermi  stores the third highest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrerrg  stores the fourth highest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrerpk  stores the fifth highest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputreltb  stores the fifth lowest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrelix  stores the fourth lowest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrelmi  stores the third lowest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrelrg  stores the second lowest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputrelpk  stores the lowest doubles of the real parts
 *               of the coefficients of the input series;
 *   inputimrtb  stores the highest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimrix  stores the second highest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimrmi  stores the third highest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimrrg  stores the fourth highest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimrpk  stores the fifth highest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimltb  stores the fifth lowest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimlix  stores the fourth lowest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimlmi  stores the third lowest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimlrg  stores the second lowest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   inputimlpk  stores the lowest doubles of the imaginary parts
 *               of the coefficients of the input series;
 *   outputrertb has space allocated for dim+1 series of degree deg,
 *               for the real highest doubles of the output;
 *   outputrerix has space allocated for dim+1 series of degree deg,
 *               for the real second highest doubles of the output;
 *   outputrermi has space allocated for dim+1 series of degree deg,
 *               for the real third highest doubles of the output;
 *   outputrerrg has space allocated for dim+1 series of degree deg,
 *               for the real fourth highest doubles of the output;
 *   outputrerpk has space allocated for dim+1 series of degree deg,
 *               for the real fifth highest doubles of the output;
 *   outputreltb has space allocated for dim+1 series of degree deg,
 *               for the real fifth lowest doubles of the output;
 *   outputrelix has space allocated for dim+1 series of degree deg,
 *               for the real fourth lowest doubles of the output;
 *   outputrelmi has space allocated for dim+1 series of degree deg,
 *               for the real third lowest doubles of the output;
 *   outputrelrg has space allocated for dim+1 series of degree deg,
 *               for the real second lowest doubles of the output;
 *   outputrelpk has space allocated for dim+1 series of degree deg,
 *               for the real lowest doubles of the output;
 *   outputimrtb has space allocated for dim+1 series of degree deg,
 *               for the imaginary highest doubles of the output;
 *   outputimrix has space allocated for dim+1 series of degree deg,
 *               for the imaginary second highest doubles of the output;
 *   outputimrmi has space allocated for dim+1 series of degree deg,
 *               for the imaginary third highest doubles of the output;
 *   outputimrrg has space allocated for dim+1 series of degree deg,
 *               for the imaginary fourth highest doubles of the output;
 *   outputimrpk has space allocated for dim+1 series of degree deg,
 *               for the imaginary fifth highest doubles of the output;
 *   outputimltb has space allocated for dim+1 series of degree deg,
 *               for the imaginary fifth lowest doubles of the output;
 *   outputimlix has space allocated for dim+1 series of degree deg,
 *               for the imaginary fourth lowest doubles of the output;
 *   outputimlmi has space allocated for dim+1 series of degree deg,
 *               for the imaginary third lowest doubles of the output;
 *   outputimlrg has space allocated for dim+1 series of degree deg,
 *               for the imaginary second lowest doubles of the output;
 *   outputimlpk has space allocated for dim+1 series of degree deg,
 *               for the imaginary lowest doubles of the output.
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
 *               derivatives and the value,
 *   outputimrpk stores the fifth highest doubles of the imaginary parts
 *               derivatives and the value,
 *   outputimltb stores the fifth lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlix stores the fourth lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlmi stores the third lowest doubles of the imaginary parts
 *               of the derivatives and the value,
 *   outputimlrg stores the second lowest doubles of the imaginary parts
 *               derivatives and the value,
 *   outputimlpk stores the lowest doubles of the imaginary parts
 *               derivatives and the value,
 *               output[idx[k]], for k from 0 to nvr, contains the
 *               derivative with respect to the variable idx[k];
 *               output[dim] contains the value of the product. */

#endif
