/* The file dbl5_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in penta double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl5_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl5_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl5_polynomials_host_h__
#define __dbl5_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl5_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double **forwardtb, double **forwardix, double **forwardmi,
   double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk,
   double **crosstb, double **crossix, double **crossmi,
   double **crossrg, double **crosspk, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in penta double precision.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   inputtb    has the highest parts of the power series
 *              for all variables in the polynomial;
 *   inputix    has the second highest parts of the power series
 *              for all variables in the polynomial;
 *   inputmi    has the middle parts of the power series
 *              for all variables in the polynomial;
 *   inputrg    has the second lowest parts of the power series
 *              for all variables in the polynomial;
 *   inputpk    has the lowest parts of the power series
 *              for all variables in the polynomial;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
 *   forwardtb  is work space for the highest doubles of nvr forward 
 *              products, forwardtb[k] can store deg+1 doubles;
 *   forwardix  is work space for the second highest doubles of nvr
 *              forward  products, forwardix[k] can store deg+1 doubles;
 *   forwardmi  is work space for the middle doubles of nvr
 *              forward products, forwardmi[k] can store deg+1 doubles;
 *   forwardrg  is work space for the second lowest doubles of nvr
 *              forward products, forwardrg[k] can store deg+1 doubles;
 *   forwardpk  is work space for the lowest doubles of nvr
 *              forward products, forwardpk[k] can store deg+1 doubles;
 *   backwardtb is work space for the highest doubles of nvr-1 backward
 *              products, backwardtb[k] can store deg+1 doubles;
 *   backwardix is work space for the second highest doubles of nvr-1
 *              backward products, backwardix[k] can store deg+1 doubles;
 *   backwardmi is work space for the middle doubles of nvr-1
 *              backward products, backwardmi[k] can store deg+1 doubles;
 *   backwardrg is work space for the second lowest doubles of nvr-1
 *              backward products, backwardrg[k] can store deg+1 doubles;
 *   backwardpk is work space for the lowest doubles of nvr-1 backward
 *              products, backwardpk[k] can store deg+1 doubles;
 *   crosstb    is work space for the highest doubles of nvr-2 cross
 *              products, crosstb[k] can store deg+1 doubles;
 *   crossix    is work space for the second highest doubles of nvr-2
 *              cross products, crossix[k] can store deg+1 doubles;
 *   crossmi    is work space for the second highest doubles of nvr-2
 *              cross products, crossmi[k] can store deg+1 doubles;
 *   crossrg    is work space for the second lowest doubles of nvr-2
 *              cross products, crossrg[k] can store for deg+1 doubles.
 *   crosspk    is work space for the lowest doubles of nvr-2 cross
 *              products, crosspk[k] can store for deg+1 doubles.
 *   verbose    if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial. */

void CPU_cmplx5_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffretb, double **cffreix, double **cffremi, 
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk,
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
   double **crossimrg, double **crossimpk, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for complex coefficients in quad double precision.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cffretb      has the highest doubles of the real parts
 *                of the coefficients, cffretb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffreix      has the second highest doubles of the real parts
 *                of the coefficients, cffreix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffremi      has the middle doubles of the real parts
 *                of the coefficients, cffremi[k] has deg+1 middle
 *                coefficients of monomial k;
 *   cffrerg      has the second lowest doubles of the real parts
 *                of the coefficients, cffrerg[k] has deg+1 second lowest
 *                coefficients of monomial k;
 *   cffrepk      has the lowest doubles of the real parts
 *                of the coefficients, cffrepk[k] has deg+1 lowest
 *                coefficients of monomial k;
 *   cffimtb      has the highest doubles of the imaginary parts
 *                of the coefficients, cffimtb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffimix      has the second highest doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimmi      has the middle doubles of the imaginary parts
 *                of the coefficients, cffimmi[k] has deg+1 middle
 *                coefficients of monomial k;
 *   cffimrg      has the second lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 second
 *                lowest coefficients of monomial k;
 *   cffimpk      has the lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 lowest
 *                coefficients of monomial k;
 *   inputretb    has the highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputreix    has the second highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputremi    has the middle doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrepk    has the second lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrepk    has the lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimtb    has the highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimix    has the second highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimmi    has the middle doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimrg    has the second lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimpk    has the lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputretb   has space for the highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputreix   has space for the second highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputremi   has space for the middle doubles of the real parts
 *                of the value and all derivatives;
 *   outputrerg   has space for the second lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrepk   has space for the lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputimtb   has space for the highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimix   has space for the second highest doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimmi   has space for the middle doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimrg   has space for the second lowest doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimpk   has space for the lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   forwardretb  has space for the highest doubles of the real parts
 *                of all nvr forward products,
 *                forwardretb[k] has space for deg+1 doubles;
 *   forwardreix  has space for the second highest doubles of the real parts
 *                of all nvr forward products,
 *                forwardreix[k] has space for deg+1 doubles;
 *   forwardremi  has space for the middle doubles of the real parts
 *                of all nvr forward products,
 *                forwardremi[k] has space for deg+1 doubles;
 *   forwardrerg  has space for the second lowest doubles of the real parts
 *                of all nvr forward products,
 *                forwardrerg[k] has space for deg+1 doubles;
 *   forwardrepk  has space for the lowest doubles of the real parts
 *                of all nvr forward products,
 *                forwardrepk[k] has space for deg+1 doubles;
 *   forwardimtb  has space for the highest doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimtb[k] has space for deg+1 doubles;
 *   forwardimix  has space for the second highest doubles of 
 *                the imaginary parts of all nvr forward products,
 *                forwardimix[k] has space for deg+1 doubles;
 *   forwardimmi  has space for the middle doubles of 
 *                the imaginary parts of all nvr forward products,
 *                forwardimmi[k] has space for deg+1 doubles;
 *   forwardimrg  has space for the second lowest doubles of the
 *                imaginary parts of all nvr forward products,
 *                forwardimrg[k] has space for deg+1 doubles;
 *   forwardimpk  has space for the lowest doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimpk[k] has space for deg+1 doubles;
 *   backwardretb has space for the highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardretb[k] has space for deg+1 doubles;
 *   backwardreix has space for the second highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardreix[k] has space for deg+1 doubles;
 *   backwardremi has space for the middle doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardremi[k] has space for deg+1 doubles;
 *   backwardrerg has space for the second lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrerg[k] has space for deg+1 doubles;
 *   backwardrepk has space for the lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrepk[k] has space for deg+1 doubles;
 *   backwardimtb has space for the highest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimtb[k] has space for deg+1 doubles;
 *   backwardimix has space for the second highest doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimix[k] has space for deg+1 doubles;
 *   backwardimmi has space for the middle doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimmi[k] has space for deg+1 doubles;
 *   backwardimrg has space for the second lowest doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimrg[k] has space for deg+1 doubles;
 *   backwardimpk has space for the lowest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimpk[k] has space for deg+1 doubles;
 *   crossretb    has space for the highest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossretb[k] has space for deg+1 doubles;
 *   crossreix    has space for the second highest doubles
 *                of the real parts of all nvr-2 cross products;
 *                crossreix[k] has space for deg+1 doubles;
 *   crossremi    has space for the middle doubles
 *                of the real parts of all nvr-2 cross products;
 *                crossremi[k] has space for deg+1 doubles;
 *   crossrerg    has space for the second lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossimrg[k] has space for deg+1 doubles;
 *   crossrepk    has space for the lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossimpk[k] has space for deg+1 doubles;
 *   crossimtb    has space for the highest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimtb[k] has space for deg+1 doubles;
 *   crossimix    has space for the second highest doubles
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimix[k] has space for deg+1 doubles;
 *   crossimmi    has space for the middle doubles
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimmi[k] has space for deg+1 doubles;
 *   crossimrg    has space for the second lowest doubles 
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimrg[k] has space for deg+1 doubles;
 *   crossimpk    has space for the lowest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimpk[k] has space for deg+1 doubles;
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputretb   has the highest doubles of the real parts,
 *   outputreix   has the second highest doubles of the real parts,
 *   outputremi   has the middle doubles of the real parts,
 *   outputrerg   has the second lowest doubles of the real parts,
 *   outputrepk   has the lowest doubles of the real parts,
 *   outputimtb   has the highest doubles of the imaginary parts,
 *   outputimix   has the second highest doubles of the imaginary parts,
 *   outputimmi   has the middle doubles of the imaginary parts,
 *   outputimrg   has the second lowest doubles of the imaginary parts,
 *   outputimpk   has the lowest doubles of the imaginary parts
 *                of derivatives and the value,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial. */

void CPU_dbl5_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second highest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   inputtb    has the highest parts of the power series
 *              for all variables in the polynomial;
 *   inputix    has the second highest parts of the power series
 *              for all variables in the polynomial;
 *   inputmi    has the middle parts of the power series
 *              for all variables in the polynomial;
 *   inputrg    has the second lowest parts of the power series
 *              for all variables in the polynomial;
 *   inputpk    has the lowest parts of the power series
 *              for all variables in the polynomial;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
 *   verbose    if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial;
 *   elapsedsec is the elapsed time in seconds. */

void CPU_cmplx5_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk, 
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk, 
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk,
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *   Evaluates and differentiates the polynomial.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cstretb      highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstreix      second highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstremx      middle deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrerg      second lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstrepk      lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstimtb      highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimix      second highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimmi      middle deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimrg      second lowest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimpk      lowest deg+1 doubles for the imaginary parts
 *                of the constant coefficient series;
 *   cffretb      has the highest doubles of the real parts
 *                of the coefficients, cffretb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffreix      has the second highest doubles of the real parts
 *                of the coefficients, cffreix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffremi      has the middle doubles of the real parts
 *                of the coefficients, cffremi[k] has deg+1 middle
 *                coefficients of monomial k;
 *   cffrerg      has the second lowest doubles of the real parts
 *                of the coefficients, cffrerg[k] has deg+1 second lowest
 *                coefficients of monomial k;
 *   cffrepk      has the lowest doubles of the real parts
 *                of the coefficients, cffrepk[k] has deg+1 lowest
 *                coefficients of monomial k;
 *   cffimtb      has the highest doubles of the imaginary parts
 *                of the coefficients, cffimtb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffimix      has the second highest doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimmi      has the middle doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 middle
 *                coefficients of monomial k;
 *   cffimrg      has the second lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 second
 *                lowest coefficients of monomial k;
 *   cffimpk      has the lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 lowest
 *                coefficients of monomial k;
 *   inputretb    has the highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputreix    has the second highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputremi    has the middle doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrepk    has the second lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrepk    has the lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimtb    has the highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimix    has the second highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimmi    has the middle doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimrg    has the second lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimpk    has the lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputretb   has space for the highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputreix   has space for the second highest doubles of the real parts
 *                of the value and all derivatives;
 *   outputremi   has space for the middle doubles of the real parts
 *                of the value and all derivatives;
 *   outputrerg   has space for the second lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputrepk   has space for the lowest doubles of the real parts
 *                of the value and all derivatives;
 *   outputimtb   has space for the highest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimix   has space for the second highest doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimmi   has space for the middle doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimrg   has space for the second lowest doubles of
 *                the imaginary parts of the value and all derivatives;
 *   outputimpk   has space for the lowest doubles of the imaginary parts
 *                of the value and all derivatives;
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputretb   has the highest doubles of the real parts,
 *   outputreix   has the second highest doubles of the real parts,
 *   outputremi   has the middle doubles of the real parts,
 *   outputrerg   has the second lowest doubles of the real parts,
 *   outputrepk   has the lowest doubles of the real parts,
 *   outputimtb   has the highest doubles of the imaginary parts,
 *   outputimix   has the second highest doubles of the imaginary parts,
 *   outputimmi   has the middle doubles of the imaginary parts,
 *   outputimrg   has the second lowest doubles of the imaginary parts,
 *   outputimpk   has the lowest doubles of the imaginary parts
 *                of derivatives and the value,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial. */

void CPU_dbl5_conv_job
 ( int deg, int nvr, int *idx,
   double *cfftb, double *cffix, double *cffmi,
   double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **forwardtb, double **forwardix, double **forwardmi,
   double **forwardrg, double **forwardpk,
   double **backwardtb, double **backwardix, double **backwardmi,
   double **backwardrg, double **backwardpk,
   double **crosstb, double **crossix, double **crossmi,
   double **crossrg, double **crosspk,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   nvr        number of variables in the monomial;
 *   idx        indices to the variables in the monomial;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   inputtb    has the highest parts of the power series
 *              for all variables in the polynomial;
 *   inputix    has the second highest parts of the power series
 *              for all variables in the polynomial;
 *   inputmi    has the middle parts of the power series
 *              for all variables in the polynomial;
 *   inputrg    has the second lowest parts of the power series
 *              for all variables in the polynomial;
 *   inputpk    has the lowest parts of the power series
 *              for all variables in the polynomial;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
 *   forwardtb  is work space for the highest doubles of nvr forward 
 *              products, forwardtb[k] can store deg+1 doubles;
 *   forwardix  is work space for the second highest doubles of nvr
 *              forward  products, forwardix[k] can store deg+1 doubles;
 *   forwardmi  is work space for the middle doubles of nvr
 *              forward  products, forwardmi[k] can store deg+1 doubles;
 *   forwardrg  iss work space for the second lowest doubles of nvr
 *              forward products, forwardrg[k] can store deg+1 doubles;
 *   forwardpk  is work space for the lowest doubles of nvr
 *              forward products, forwardpk[k] can store deg+1 doubles;
 *   backwardtb is work space for the highest doubles of nvr-1 backward
 *              products, backwardtb[k] can store deg+1 doubles;
 *   backwardix is work space for the second highest doubles of nvr-1
 *              backward products, backwardix[k] can store deg+1 doubles;
 *   backwardmi is work space for the middle doubles of nvr-1
 *              backward products, backwardmi[k] can store deg+1 doubles;
 *   backwardrg is work space for the second lowest doubles of nvr-1
 *              backward products, backwardrg[k] can store deg+1 doubles;
 *   backwardpk is work space for the lowest doubles of nvr-1 backward
 *              products, backwardpk[k] can store deg+1 doubles;
 *   crosstb    is work space for the highest doubles of nvr-2 cross
 *              products, crosstb[k] can store deg+1 doubles;
 *   crossix    is work space for the second highest doubles of nvr-2
 *              cross products, crossix[k] can store deg+1 doubles;
 *   crossmi    is work space for the middle doubles of nvr-2
 *              cross products, crossmi[k] can store deg+1 doubles;
 *   crossrg    is work space for the second lowest doubles of nvr-2
 *              cross products, crossrg[k] can store for deg+1 doubles;
 *   crosspk    is work space for the lowest doubles of nvr-2 cross
 *              products, crosspk[k] can store for deg+1 doubles;
 *   job        defines the convolution job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardtb  are the updated highest parts of forward products;
 *   forwardix  are the updated second highest parts of forward products;
 *   forwardmi  are the updated middle parts of forward products;
 *   forwardrg  are the updated second lowest parts of forward products;
 *   forwardpk  are the updated lowest parts forward products;
 *   backwardtb are the updated highest parts of backward products;
 *   backwardix are the updated second highest parts of backward products;
 *   backwardmi are the updated middle parts of backward products;
 *   backwardrg are the updated second lowest parts of backward products;
 *   backwardpk are the updated lowest parts backward products;
 *   crosstb    are the updated highest parts of cross products;
 *   crossix    are the updated second highest parts of cross products;
 *   crossmi    are the updated middle parts of cross products;
 *   crossrg    are the updated second lowest parts of cross products;
 *   crosspk    are the updated lowest parts cross products. */

void CPU_cmplx5_conv_job
 ( int deg, int nvr, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
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
   double **crossimrg, double **crossimpk,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   nvr          number of variables in the monomial;
 *   idx          indices to the variables in the monomial;
 *   cffretb      has the highest doubles of the real parts
 *                of the coefficients, cffretb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffreix      has the second highest doubles of the real parts
 *                of the coefficients, cffreix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffremi      has the middle doubles of the real parts
 *                of the coefficients, cffremi[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffrerg      has the second lowest doubles of the real parts
 *                of the coefficients, cffrerg[k] has deg+1 second lowest
 *                coefficients of monomial k;
 *   cffrepk      has the lowest doubles of the real parts
 *                of the coefficients, cffrepk[k] has deg+1 lowest
 *                coefficients of monomial k;
 *   cffimtb      has the highest doubles of the imaginary parts
 *                of the coefficients, cffimtb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffimix      has the second highest doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimmi      has the middle doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimrg      has the second lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 second
 *                lowest coefficients of monomial k;
 *   cffimpk      has the lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 lowest
 *                coefficients of monomial k;
 *   inputretb    has the highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputreix    has the second highest doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputremi    has the middle doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrerg    has the second lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrepk    has the lowest doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimtb    has the highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimix    has the second highest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimmi    has the middle doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimrg    has the second lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimpk    has the lowest doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   forwardretb  has space for the highest doubles of the real parts
 *                of all nvr forward products,
 *                forwardretb[k] has space for deg+1 doubles;
 *   forwardreix  has space for the second highest doubles of the real parts
 *                of all nvr forward products,
 *                forwardreix[k] has space for deg+1 doubles;
 *   forwardremi  has space for the second highest doubles of the real parts
 *                of all nvr forward products,
 *                forwardremi[k] has space for deg+1 doubles;
 *   forwardrerg  has space for the second lowest doubles of the real parts
 *                of all nvr forward products,
 *                forwardrerg[k] has space for deg+1 doubles;
 *   forwardrepk  has space for the lowest doubles of the real parts
 *                of all nvr forward products,
 *                forwardrepk[k] has space for deg+1 doubles;
 *   forwardimtb  has space for the highest doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimtb[k] has space for deg+1 doubles;
 *   forwardimix  has space for the second highest doubles of 
 *                the imaginary parts of all nvr forward products,
 *                forwardimix[k] has space for deg+1 doubles;
 *   forwardimmi  has space for the second highest doubles of 
 *                the imaginary parts of all nvr forward products,
 *                forwardimmi[k] has space for deg+1 doubles;
 *   forwardimrg  has space for the second lowest doubles of the
 *                imaginary parts of all nvr forward products,
 *                forwardimrg[k] has space for deg+1 doubles;
 *   forwardimpk  has space for the lowest doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimpk[k] has space for deg+1 doubles;
 *   backwardretb has space for the highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardretb[k] has space for deg+1 doubles;
 *   backwardreix has space for the second highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardreix[k] has space for deg+1 doubles;
 *   backwardremi has space for the middle doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardremi[k] has space for deg+1 doubles;
 *   backwardrerg has space for the second lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrerg[k] has space for deg+1 doubles;
 *   backwardrepk has space for the lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrepk[k] has space for deg+1 doubles;
 *   backwardimtb has space for the highest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimtb[k] has space for deg+1 doubles;
 *   backwardimix has space for the second highest doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimix[k] has space for deg+1 doubles;
 *   backwardimmi has space for the middle doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimmi[k] has space for deg+1 doubles;
 *   backwardimrg has space for the second lowest doubles
 *                of the imaginary parts of all nvr-2 backward products;
 *                backwardimrg[k] has space for deg+1 doubles;
 *   backwardimpk has space for the lowest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimpk[k] has space for deg+1 doubles;
 *   crossretb    has space for the highest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossretb[k] has space for deg+1 doubles;
 *   crossreix    has space for the second highest doubles
 *                of the real parts of all nvr-2 cross products;
 *                crossreix[k] has space for deg+1 doubles;
 *   crossremi    has space for the middle doubles
 *                of the real parts of all nvr-2 cross products;
 *                crossremi[k] has space for deg+1 doubles;
 *   crossrerg    has space for the second lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossimrg[k] has space for deg+1 doubles;
 *   crossrepk    has space for the lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossimpk[k] has space for deg+1 doubles;
 *   crossimtb    has space for the highest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimtb[k] has space for deg+1 doubles;
 *   crossimix    has space for the second highest doubles
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimix[k] has space for deg+1 doubles;
 *   crossimmi    has space for the middle doubles
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimmi[k] has space for deg+1 doubles;
 *   crossimrg    has space for the second lowest doubles 
 *                of the imaginary parts of all nvr-2 cross products;
 *                crossimrg[k] has space for deg+1 doubles;
 *   crossimpk    has space for the lowest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimpk[k] has space for deg+1 doubles;
 *   job          defines the convolution job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardretb  are the updated highest doubles of the real parts
 *                of the forward products;
 *   forwardreix  are the updated second highest doubles of the real parts
 *                of the forward products;
 *   forwardremi  are the updated middle doubles of the real parts
 *                of the forward products;
 *   forwardrerg  are the updated second lowest doubles of the real parts
 *                of the forward products;
 *   forwardrepk  are the updated lowest doubles of the real parts
 *                of the forward products;
 *   forwardimtb  are the updated highest doubles of the imaginary parts
 *                of the forward products;
 *   forwardimix  are the updated second highest doubles of the
 *                imaginary parts of the forward products;
 *   forwardimmi  are the updated middle doubles of the
 *                imaginary parts of the forward products;
 *   forwardimrg  are the updated second lowest doubles of the imaginary parts
 *                of the forward products;
 *   forwardimpk  are the updated lowest doubles of the imaginary parts
 *                of the forward products;
 *   backwardretb are the updated highest doubles of the real parts 
 *                of the backward products;
 *   backwardreix are the updated second highest doubles of the real parts 
 *                of the backward products;
 *   backwardremi are the updated middle doubles of the real parts 
 *                of the backward products;
 *   backwardrerg are the updated second lowest doubles of the real parts 
 *                of the backward products;
 *   backwardrepk are the updated lowest doubles of the real parts 
 *                of the backward products;
 *   backwardimtb are the updated highest doubles of the imaginary parts 
 *                of the backward products;
 *   backwardimix are the updated second highest doubles of
 *                the imaginary parts of the backward products;
 *   backwardimmi are the updated middle doubles of
 *                the imaginary parts of the backward products;
 *   backwardimrg are the updated second lowest doubles of
 *                the imaginary parts of the backward products;
 *   backwardimpk are the updated lowest doubles of the imaginary parts 
 *                of the backward products;
 *   crossretb    are the updated highest doubles of the real parts
 *                of the cross products;
 *   crossreix    are the updated second highest doubles of the real parts
 *                of the cross products;
 *   crossremi    are the updated middle doubles of the real parts
 *                of the cross products;
 *   crossrerg    are the updated second lowest doubles of the real parts
 *                of the cross products;
 *   crossrepk    are the updated lowest doubles of the real parts
 *                of the cross products;
 *   crossimtb    are the updated highest doubles of the imaginary parts
 *                of the cross products;
 *   crossimix    are the updated second highest doubles of
 *                the imaginary parts of the cross products;
 *   crossimmi    are the updated middle doubles of
 *                the imaginary parts of the cross products;
 *   crossimrg    are the updated second lowest doubles of
 *                the imaginary parts of the cross products;
 *   crossimpk    are the updated lowest doubles of the imaginary parts
 *                of the cross products. */

void CPU_dbl5_add_job
 ( int deg,
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second highest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   forwardtb  are all highest parts of computed forward products,
 *   forwardix  are all second highest parts of computed forward products,
 *   forwardmi  are all middle parts of computed forward products,
 *   forwardrg  are all second lowest parts of computed forward products,
 *   forwardpk  are all lowest parts of computed forward products,
 *   backwardtb are all highest parts of computed backward products;
 *   backwardix are all second highest parts of computed backward products;
 *   backwardmi are all middle parts of computed backward products;
 *   backwardrg are all second lowest parts of computed backward products;
 *   backwardpk are all lowest parts of computed backward products;
 *   crosstb    are all highest parts of computed cross products;
 *   crossix    are all second highest parts of computed cross products;
 *   crossmi    are all middle parts of computed cross products;
 *   crossrg    are all second lowest parts of computed cross products;
 *   crosspk    are all lowest parts of computed cross products;
 *   job        defines the addition job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardtb  are the updated highest parts of forward products;
 *   forwardix  are the updated second highest parts of forward products;
 *   forwardmi  are the updated middle parts of forward products;
 *   forwardrg  are the updated second lowest parts of forward products;
 *   forwardpk  are the updated lowest parts forward products;
 *   backwardtb are the updated highest parts of backward products;
 *   backwardix are the updated second highest parts of backward products;
 *   backwardmi are the updated middle parts of backward products;
 *   backwardrg are the updated second lowest parts of backward products;
 *   backwardpk are the updated lowest parts backward products;
 *   crosstb    are the updated highest parts of cross products;
 *   crossix    are the updated second highest parts of cross products;
 *   crossmi    are the updated middle parts of cross products;
 *   crossrg    are the updated second lowest parts of cross products;
 *   crosspk    are the updated lowest parts cross products. */

void CPU_cmplx5_add_job
 ( int deg, double *cstretb, double *cstreix, double *cstremi,
            double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk,
   double ***forwardretb, double ***forwardreix, double ***forwardremi,
   double ***forwardrerg, double ***forwardrepk,
   double ***forwardimtb, double ***forwardimix, double ***forwardimmi,
   double ***forwardimrg, double ***forwardimpk,
   double ***backwardretb, double ***backwardreix, double ***backwardremi,
   double ***backwardrerg, double ***backwardrepk, 
   double ***backwardimtb, double ***backwardimix, double ***backwardimmi,
   double ***backwardimrg, double ***backwardimpk, 
   double ***crossretb, double ***crossreix, double ***crossremi,
   double ***crossrerg, double ***crossrepk,
   double ***crossimtb, double ***crossimix, double ***crossimmi,
   double ***crossimrg, double ***crossimpk,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   cstretb      highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstreix      second highest deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstremi      middle deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrerg      second lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstrepk      lowest deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstimtb      highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimix      second highest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimmi      middle deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimrg      second lowest deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimpk      lowest deg+1 doubles for the imaginary parts
 *                of the constant coefficient series;
 *   cffretb      has the highest doubles of the real parts
 *                of the coefficients, cffretb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffreix      has the second highest doubles of the real parts
 *                of the coefficients, cffreix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffremi      has the middle doubles of the real parts
 *                of the coefficients, cffremi[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffrerg      has the second lowest doubles of the real parts
 *                of the coefficients, cffrerg[k] has deg+1 second lowest
 *                coefficients of monomial k;
 *   cffrepk      has the lowest doubles of the real parts
 *                of the coefficients, cffrepk[k] has deg+1 lowest
 *                coefficients of monomial k;
 *   cffimtb      has the highest doubles of the imaginary parts
 *                of the coefficients, cffimtb[k] has deg+1 highest
 *                coefficients of monomial k;
 *   cffimix      has the second highest doubles of the imaginary parts
 *                of the coefficients, cffimix[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimmi      has the middle doubles of the imaginary parts
 *                of the coefficients, cffimmi[k] has deg+1 second highest
 *                coefficients of monomial k;
 *   cffimrg      has the second lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 second
 *                lowest coefficients of monomial k;
 *   cffimpk      has the lowest doubles of the imaginary parts
 *                of the coefficient, cffimpk[k] has the deg+1 lowest
 *                coefficients of monomial k;
 *   forwardretb  computed highest doubles of the real parts
 *                of all nvr forward products;
 *   forwardreix  computed second highest doubles of the real parts
 *                of all nvr forward products;
 *   forwardremi  computed middle doubles of the real parts
 *                of all nvr forward products;
 *   forwardrerg  computed second lowest doubles of the real parts
 *                of all nvr forward products;
 *   forwardrepk  computed lowest doubles of the real parts
 *                of all nvr forward products;
 *   forwardimtb  computed highest doubles of the imaginary parts
 *                of all nvr forward products;
 *   forwardimix  computed second highest doubles of the imaginary parts
 *                of all nvr forward products;
 *   forwardimmi  computed middle doubles of the imaginary parts
 *                of all nvr forward products;
 *   forwardimrg  computed second lowest doubles of the imaginary parts
 *                of all nvr forward products;
 *   forwardimpk  computed lowest doubles of the imaginary parts
 *                of all nvr forward products;
 *   backwardretb computed highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *   backwardreix computed second highest doubles of the real parts 
 *                of all nvr-2 backward products;
 *   backwardremi computed middle doubles of the real parts 
 *                of all nvr-2 backward products;
 *   backwardrerg computed second lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *   backwardrepk computed lowest doubles of the real parts 
 *                of all nvr-2 backward products;
 *   backwardimtb computed highest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *   backwardimix computed second highest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *   backwardimmi computed middle doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *   backwardimrg computed second lowest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *   backwardimpk computed lowest doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *   crossretb    computed highest doubles of the real parts
 *                of all nvr-2 cross products;
 *   crossreix    computed second highest doubles of the real parts
 *                of all nvr-2 cross products;
 *   crossremi    computed middle doubles of the real parts
 *                of all nvr-2 cross products;
 *   crossrerg    computed second lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *   crossrepk    computed lowest doubles of the real parts
 *                of all nvr-2 cross products;
 *   crossimtb    computed highest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *   crossimix    computed second highest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *   crossimmi    computed middle doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *   crossimrg    computed second lowest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *   crossimpk    computed lowest doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *   job          defines the addition job;
 *   verbose      if true, then is verbose.

 * ON RETURN :
 *   forwardretb  are the updated highest doubles of the real parts
 *                of the forward products;
 *   forwardreix  are the updated second highest doubles of the real parts
 *                of the forward products;
 *   forwardremi  are the updated middle doubles of the real parts
 *                of the forward products;
 *   forwardrerg  are the updated second lowest doubles of the real parts
 *                of the forward products;
 *   forwardrepk  are the updated lowest doubles of the real parts
 *                of the forward products;
 *   forwardimtb  are the updated highest doubles of the imaginary parts
 *                of the forward products;
 *   forwardimix  are the updated second highest doubles of the
 *                imaginary parts of the forward products;
 *   forwardimmi  are the updated middle doubles of the
 *                imaginary parts of the forward products;
 *   forwardimrg  are the updated second lowest doubles of the
 *                imaginary parts of the forward products;
 *   forwardimpk  are the updated lowest doubles of the imaginary parts
 *                of the forward products;
 *   backwardretb are the updated highest doubles of the real parts 
 *                of the backward products;
 *   backwardreix are the updated second highest doubles of the real parts 
 *                of the backward products;
 *   backwardremi are the updated middle doubles of the real parts 
 *                of the backward products;
 *   backwardrerg are the updated second lowest doubles of the real parts 
 *                of the backward products;
 *   backwardrepk are the updated lowest doubles of the real parts 
 *                of the backward products;
 *   backwardimtb are the updated highest doubles of the imaginary parts 
 *                of the backward products;
 *   backwardimix are the updated second highest doubles of the
 *                imaginary parts of the backward products;
 *   backwardimmi are the updated middle doubles of the
 *                imaginary parts of the backward products;
 *   backwardimrg are the updated second lowest doubles of the
 *                imaginary parts of the backward products;
 *   backwardimpk are the updated lowest doubles of the imaginary parts 
 *                of the backward products;
 *   crossretb    are the updated highest doubles of the real parts
 *                of the cross products;
 *   crossreix    are the updated second highest doubles of the real parts
 *                of the cross products;
 *   crossremi    are the updated middle doubles of the real parts
 *                of the cross products;
 *   crossrerg    are the updated second lowest doubles of the real parts
 *                of the cross products;
 *   crossrepk    are the updated lowest doubles of the real parts
 *                of the cross products;
 *   crossimtb    are the updated highest doubles of the
 *                imaginary parts of the cross products;
 *   crossimix    are the updated second highest doubles of the
 *                imaginary parts of the cross products;
 *   crossimmi    are the updated middle doubles of the
 *                imaginary parts of the cross products;
 *   crossimrg    are the updated lowest doubles of the imaginary parts
 *                of the cross products;
 *   crossimpk    are the updated lowest doubles of the imaginary parts
 *                of the cross products. */

void CPU_dbl5_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second highest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   forwardtb  are all highest parts of computed forward products;
 *   forwardix  are all second highest parts of computed forward products;
 *   forwardmi  are all middle parts of computed forward products;
 *   forwardrg  are all second lowest parts of computed forward products;
 *   forwardpk  are all lowest parts of computed forward products;
 *   backwardtb are all highest parts of computed backward products;
 *   backwardix are all second highest parts of computed backward products;
 *   backwardmi are all middle parts of computed backward products;
 *   backwardrg are all second lowest parts of computed backward products;
 *   backwardpk are all lowest parts of computed backward products;
 *   crosstb    are all highest parts of computed cross products;
 *   crossix    are all second highest parts of computed cross products;
 *   crossmi    are all middle parts of computed cross products;
 *   crossrg    are all second lowest parts of computed cross products;
 *   crosspk    are all lowest parts of computed cross products.
 *
 * ON RETURN :
 *   outputtb   highest parts of the values and all derivatives;
 *   outputix   second highest parts of the values and all derivatives;
 *   outputmi   middle parts of the values and all derivatives;
 *   outputrg   second lowest parts of the values and all derivatives;
 *   outputpk   lowest parts of the values and all derivatives. */

void CPU_dbl5_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   double ***forwardtb, double ***forwardix, double ***forwardmi,
   double ***forwardrg, double ***forwardpk,
   double ***backwardtb, double ***backwardix, double ***backwardmi,
   double ***backwardrg, double ***backwardpk, 
   double ***crosstb, double ***crossix, double ***crossmi,
   double ***crossrg, double ***crosspk,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second highest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   forwardtb  are all highest parts of computed forward products;
 *   forwardix  are all second highest parts of computed forward products;
 *   forwardmi  are all middle parts of computed forward products;
 *   forwardrg  are all second lowest parts of computed forward products;
 *   forwardpk  are all lowest parts of computed forward products;
 *   backwardtb are all highest parts of computed backward products;
 *   backwardix are all second highest parts of computed backward products;
 *   backwardmx are all middle parts of computed backward products;
 *   backwardrg are all second lowest parts of computed backward products;
 *   backwardpk are all lowest parts of computed backward products;
 *   crosstb    are all highest parts of computed cross products;
 *   crossix    are all second highest parts of computed cross products;
 *   crossmi    are all middle parts of computed cross products;
 *   crossrg    are all second lowest parts of computed cross products;
 *   crosspk    are all lowest parts of computed cross products;
 *   jobs       defines the addition jobs;
 *   verbose    if true, then output is written.
 *
 * ON RETURN :
 *   outputtb   highest parts of the values and all derivatives;
 *   outputix   second highest parts of the values and all derivatives;
 *   outputmi   middle parts of the values and all derivatives;
 *   outputrg   second lowest parts of the values and all derivatives;
 *   outputpk   lowest parts of the values and all derivatives. */

void CPU_dbl5_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk, 
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csttb      highest parts of constant coefficient series;
 *   cstix      second highest parts of constant coefficient series;
 *   cstmi      middle parts of constant coefficient series;
 *   cstrg      second lowest parts of constant coefficient series;
 *   cstpk      lowest parts of constant coefficient series;
 *   cfftb      cfftb[k] has deg+1 doubles for the highest parts
 *              of the coefficient series of monomial k;
 *   cffix      cffix[k] has deg+1 doubles for the second highest parts
 *              of the coefficient series of monomial k;
 *   cffmi      cffmi[k] has deg+1 doubles for the middle parts
 *              of the coefficient series of monomial k;
 *   cffrg      cffrg[k] has deg+1 doubles for the second lowest parts
 *              of the coefficient series of monomial k;
 *   cffpk      cffpk[k] has deg+1 doubles for the lowest parts
 *              of the coefficient series of monomial k;
 *   inputtb    has the highest parts of the power series
 *              for all variables in the polynomial;
 *   inputix    has the second highest parts of the power series
 *              for all variables in the polynomial;
 *   inputmi    has the middle parts of the power series
 *              for all variables in the polynomial;
 *   inputrg    has the second lowest parts of the power series
 *              for all variables in the polynomial;
 *   inputpk    has the lowest parts of the power series
 *              for all variables in the polynomial;
 *   outputtb   has space allocated for dim+1 series of degree deg;
 *   outputix   has space allocated for dim+1 series of degree deg;
 *   outputmi   has space allocated for dim+1 series of degree deg;
 *   outputrg   has space allocated for dim+1 series of degree deg;
 *   outputpk   has space allocated for dim+1 series of degree deg;
 *   cnvjobs    convolution jobs organized in layers;
 *   addjobs    addition jobs organized in layers;
 *   verbose    if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputtb   has the highest parts of derivatives and the value,
 *              outputtb[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputtb[dim] contains the value of the polynomial;
 *   outputix   has the second highest parts of derivatives and the value,
 *              outputix[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputix[dim] contains the value of the polynomial;
 *   outputmi   has the middle parts of derivatives and the value,
 *              outputmi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputmi[dim] contains the value of the polynomial;
 *   outputrg   has the second lowest parts of derivatives and the value,
 *              outputrg[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputrg[dim] contains the value of the polynomial;
 *   outputpk   has the lowest parts of derivatives and the value,
 *              outputpk[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputpk[dim] contains the value of the polynomial;
 *   elapsedsec is the elapsed time in seconds. */

#endif
