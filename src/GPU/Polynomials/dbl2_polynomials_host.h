/* The file dbl2_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in double double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl2_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl2_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl2_polynomials_host_h__
#define __dbl2_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl2_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   double **forwardhi, double **forwardlo,
   double **backwardhi, double **backwardlo,
   double **crosshi, double **crosslo, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in double double precision.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   cffhi      cffhi[k] has deg+1 doubles for the high parts of the
 *              coefficient series of monomial k;
 *   cfflo      cfflo[k] has deg+1 doubles for the low parts of the
 *              coefficient series of monomial k;
 *   inputhi    has the high parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   inputlo    has the low parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   outputhi   has space for high parts of the value and all derivatives;
 *   outputlo   has space for low parts of the value and all derivatives;
 *   forwardhi  has space for high parts of all nvr forward products,
 *              forwardhi[k] has space for deg+1 doubles;
 *   forwardlo  has space for low parts of all nvr forward products,
 *              forwardlo[k] has space for deg+1 doubles;
 *   backwardhi has space for high parts of all nvr-2 backward products;
 *              backwardhi[k] has space for deg+1 doubles;
 *   backwardlo has space for low parts of all nvr-2 backward products;
 *              backwardlo[k] has space for deg+1 doubles;
 *   crosshi    has space for high parts of all nvr-2 cross products;
 *              crosshi[k] has space for deg+1 doubles;
 *   crosslo    has space for low parts of all nvr-2 cross products;
 *              crosslo[k] has space for deg+1 doubles;
 *   verbose    if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhi   has the high parts of derivatives and the value,
 *              outputhi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputhi[dim] contains the value of the polynomial;
 *   outputlo   has the low parts of derivatives and the value,
 *              outputlo[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlo[dim] contains the value of the polynomial. */

void CPU_cmplx2_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   double **forwardrehi, double **forwardrelo,
   double **forwardimhi, double **forwardimlo,
   double **backwardrehi, double **backwardrelo,
   double **backwardimhi, double **backwardimlo,
   double **crossrehi, double **crossrelo,
   double **crossimhi, double **crossimlo, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for complex coefficients in double double precision.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   inputrehi    has the high doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelo    has the low doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhi    has the high doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlo    has the low doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputrehi   has space for the high doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelo   has space for the low doubles of the real parts
 *                of the value and all derivatives;
 *   outputimhi   has space for the high doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlo   has space for the low doubles of the imaginary parts
 *                of the value and all derivatives;
 *   forwardrehi  has space for the high doubles of the real parts
 *                of all nvr forward products,
 *                forwardrehi[k] has space for deg+1 doubles;
 *   forwardrelo  has space for the low doubles of the real parts
 *                of all nvr forward products,
 *                forwardrelo[k] has space for deg+1 doubles;
 *   forwardimhi  has space for the high doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimhi[k] has space for deg+1 doubles;
 *   forwardimlo  has space for the low doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimlo[k] has space for deg+1 doubles;
 *   backwardrehi has space for the high doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrehi[k] has space for deg+1 doubles;
 *   backwardrelo has space for the low doubles of the real parts 
 *                of all nvr-2 backward products;
 *                backwardrelo[k] has space for deg+1 doubles;
 *   backwardimhi has space for the high doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimhi[k] has space for deg+1 doubles;
 *   backwardimlo has space for the low doubles of the imaginary parts 
 *                of all nvr-2 backward products;
 *                backwardimlo[k] has space for deg+1 doubles;
 *   crossrehi    has space for the high doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossrehi[k] has space for deg+1 doubles;
 *   crossrelo    has space for the low doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossrelo[k] has space for deg+1 doubles;
 *   crossimhi    has space for the high doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimhi[k] has space for deg+1 doubles;
 *   crossimlo    has space for the low doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimlo[k] has space for deg+1 doubles;
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputrehi   has the high doubles of the real parts,
 *   outputrelo   has the low doubles of the real parts,
 *   outputimhi   has the high doubles of the imaginary parts,
 *   outputimlo   has the low doubles of the imaginary parts
 *                of derivatives and the value,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial. */

void CPU_dbl2_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstlo,
   double **cffhi, double **cfflo,
   double **inputhi, double **inputlo, 
   double **outputhi, double **outputlo,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *   Evaluates and differentiates the polynomial.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csthi      high deg+1 doubles for the constant coefficient series;
 *   cstlo      low deg+1 doubles for the constant coefficient series;
 *   cffhi      cffhi[k] has deg+1 doubles for the high parts of the
 *              coefficient series of monomial k;
 *   cfflo      cfflo[k] has deg+1 doubles for the low parts of the
 *              coefficient series of monomial k;
 *   inputhi    has the high parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   inputlo    has the low parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   outputhi   has space for high parts of the value and all derivatives;
 *   outputlo   has space for low parts of the value and all derivatives;
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   outputhi   has the high parts of derivatives and the value,
 *              outputhi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputhi[dim] contains the value of the polynomial;
 *   outputlo   has the low parts of derivatives and the value,
 *              outputlo[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlo[dim] contains the value of the polynomial;
 *   elapsedsec is the elapsed time in seconds. */

void CPU_cmplx2_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo, 
   double **inputimhi, double **inputimlo, 
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   double *elapsedsec, int vrblvl );
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
 *   cstrehi      high deg+1 doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low deg+1 doubles for the real parts
 *                of the constant coefficient series;
 *   cstimhi      high deg+1 doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low deg+1 doubles for the imaginary parts
 *                of the constant coefficient series;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   inputrehi    has the high doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelo    has the low doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhi    has the high doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlo    has the low doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputrehi   has space for the high doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelo   has space for the low doubles of the real parts
 *                of the value and all derivatives;
 *   outputimhi   has space for the high doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlo   has space for the low doubles of the imaginary parts
 *                of the value and all derivatives;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputrehi   has the high doubles of the real parts,
 *   outputrelo   has the low doubles of the real parts,
 *   outputimhi   has the high doubles of the imaginary parts,
 *   outputimlo   has the low doubles of the imaginary parts
 *                of derivatives and the value,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

void CPU_dbl2_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhi, double *cfflo,
   double **inputhi, double **inputlo,
   double **forwardhi, double **forwardlo,
   double **backwardhi, double **backwardlo,
   double **crosshi, double **crosslo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job, on real data.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   nvr        number of variables in the monomial;
 *   idx        indices to the variables in the monomial;
 *   cffhi      high parts of the coefficient series of the monomial;
 *   cfflo      low parts of the coefficient series of the monomial;
 *   inputhi    has the high parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   inputlo    has the low parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   forwardhi  has space for high parts of all nvr forward products,
 *              forwardhi[k] has space for deg+1 doubles;
 *   forwardlo  has space for low parts of all nvr forward products,
 *              forwardlo[k] has space for deg+1 doubles;
 *   backwardhi has space for high parts of all nvr-2 backward products;
 *              backwardhi[k] has space for deg+1 doubles;
 *   backwardlo has space for low parts of all nvr-2 backward products;
 *              backwardlo[k] has space for deg+1 doubles;
 *   crosshi    has space for high parts of all nvr-2 cross products;
 *              crosshi[k] has space for deg+1 doubles;
 *   crosslo    has space for low parts of all nvr-2 cross products;
 *              crosslo[k] has space for deg+1 doubles;
 *   job        defines the convolution job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhi  are the updated high parts of forward products;
 *   forwardlo  are the updated low parts of forward products;
 *   backwardhi are the updated high parts of backward products;
 *   backwardlo are the updated low parts of backward products;
 *   crosshi    are the updated high parts of cross products;
 *   crosslo    are the updated low parts of cross products. */

void CPU_cmplx2_conv_job
 ( int deg, int nvr, int *idx,
   double *cffrehi, double *cffrelo,
   double *cffimhi, double *cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **forwardrehi, double **forwardrelo,
   double **forwardimhi, double **forwardimlo,
   double **backwardrehi, double **backwardrelo,
   double **backwardimhi, double **backwardimlo,
   double **crossrehi, double **crossrelo,
   double **crossimhi, double **crossimlo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   nvr          number of variables in the monomial;
 *   idx          indices to the variables in the monomial;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   inputrehi    has the high doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelo    has the low doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhi    has the high doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlo    has the low doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   forwardrehi  has space for the high doubles of the real parts
 *                of all nvr forward products,
 *                forwardrehi[k] has space for deg+1 doubles;
 *   forwardrelo  has space for the low doubles of the real parts
 *                of all nvr forward products,
 *                forwardrelo[k] has space for deg+1 doubles;
 *   forwardimhi  has space for the high doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimhi[k] has space for deg+1 doubles;
 *   forwardimlo  has space for the low doubles of the imaginary parts
 *                of all nvr forward products,
 *                forwardimlo[k] has space for deg+1 doubles;
 *   backwardrehi has space for the high doubles of the real parts
 *                of all nvr-2 backward products;
 *                backwardrehi[k] has space for deg+1 doubles;
 *   backwardrelo has space for the low doubles of the real parts
 *                of all nvr-2 backward products;
 *                backwardrelo[k] has space for deg+1 doubles;
 *   backwardimhi has space for the high doubles of the imaginary parts
 *                of all nvr-2 backward products;
 *                backwardimhi[k] has space for deg+1 doubles;
 *   backwardimlo has space for the low doubles of the imaginary parts
 *                of all nvr-2 backward products;
 *                backwardimlo[k] has space for deg+1 doubles;
 *   crossrehi    has space for the high doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossrehi[k] has space for deg+1 doubles;
 *   crossrelo    has space for the low doubles of the real parts
 *                of all nvr-2 cross products;
 *                crossrelo[k] has space for deg+1 doubles;
 *   crossimhi    has space for the high doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimhi[k] has space for deg+1 doubles;
 *   crossimlo    has space for the low doubles of the imaginary parts
 *                of all nvr-2 cross products;
 *                crossimlo[k] has space for deg+1 doubles;
 *   job          defines the convolution job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrehi  are the updated high doubles of the real parts
 *                of the forward products;
 *   forwardrelo  are the updated low doubles of the real parts
 *                of the forward products;
 *   forwardimhi  are the updated high doubles of the imaginary parts
 *                of the forward products;
 *   forwardimlo  are the updated low doubles of the imaginary parts
 *                of the forward products;
 *   backwardrehi are the updated high doubles of the real parts 
 *                of the backward products;
 *   backwardrelo are the updated low doubles of the real parts 
 *                of the backward products;
 *   backwardimhi are the updated high doubles of the imaginary parts 
 *                of the backward products;
 *   backwardimlo are the updated low doubles of the imaginary parts 
 *                of the backward products;
 *   crossrehi    are the updated high doubles of the real parts
 *                of the cross products;
 *   crossrelo    are the updated low doubles of the real parts
 *                of the cross products;
 *   crossimhi    are the updated high doubles of the imaginary parts
 *                of the cross products;
 *   crossimlo    are the updated low doubles of the imaginary parts
 *                of the cross products. */

void CPU_dbl2_add_job
 ( int deg, double *csthi, double *cstlo,
   double **cffhi, double **cfflo,
   double ***forwardhi, double ***forwardlo,
   double ***backwardhi, double ***backwardlo, 
   double ***crosshi, double ***crosslo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job, on real data.
 *
 * ON ENTRY :
 *   deg        degree of the series;
 *   csthi      high parts of constant coefficient series;
 *   cstlo      low parts of constant coefficient series;
 *   cffhi      high parts of coefficients of the monomials;
 *   cfflo      low parts of coefficients of the monomials;
 *   forwardhi  all computed high parts of forward products,
 *   forwardlo  all computed low parts of forward products,
 *   backwardhi all computed high parts of backward products;
 *   backwardlo all computed low parts of backward products;
 *   crosshi    all computed high parts of cross products;
 *   crosslo    all computed low parts of cross products;
 *   job        defines the addition job;
 *   verbose    if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhi  are the updated high parts of forward products;
 *   forwardlo  are the updated low parts forward products;
 *   backwardhi are the updated high parts of backward products;
 *   backwardlo are the updated low parts backward products;
 *   crosshi    are the updated high parts of cross products;
 *   crosslo    are the updated low parts cross products. */

void CPU_cmplx2_add_job
 ( int deg, double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double ***forwardrehi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardrelo, 
   double ***backwardimhi, double ***backwardimlo, 
   double ***crossrehi, double ***crossrelo,
   double ***crossimhi, double ***crossimlo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job, on complex data.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   cstrehi      high doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low doubles of the real parts
 *                of the constant coefficient series;
 *   cstimhi      high doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low doubles of the imaginary parts
 *                of the constant coefficient series;
 *   forwardrehi  are all computed high doubles of the real parts
 *                of forward products;
 *   forwardrelo  are all computed low doubles of the real parts
 *                of forward products;
 *   forwardimhi  are all computed high doubles of the imaginary parts
 *                of forward products;
 *   forwardimlo  are all computed low doubles of the imaginary parts
 *                of forward products;
 *   backwardrehi are all computed high doubles of the real parts
 *                of the backward products;
 *   backwardrelo are all computed low doubles of the real parts
 *                of the backward products;
 *   backwardimhi are all computed high doubles of the imaginary parts
 *                of the backward products;
 *   backwardimlo are all computed low doubles of the imaginary parts
 *                of the backward products;
 *   crossrehi    are all computed high doubles of the real parts
 *                of the cross products;
 *   crossrelo    are all computed low doubles of the real parts
 *                of the cross products;
 *   crossimhi    are all computed high doubles of the imaginary parts
 *                of the cross products;
 *   crossimlo    are all computed low doubles of the imaginary parts
 *                of the cross products;
 *   job          defines the addition job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardrehi  are the updated high doubles of the real parts
 *                of the forward products;
 *   forwardrelo  are the updated low doubles of the real parts
 *                of the forward products;
 *   forwardimhi  are the updated high doubles of the imaginary parts
 *                of the forward products;
 *   forwardimlo  are the updated low doubles of the imaginary parts
 *                of the forward products;
 *   backwardrehi are the updated high doubles of the real parts 
 *                of the backward products;
 *   backwardrelo are the updated low doubles of the real parts 
 *                of the backward products;
 *   backwardimhi are the updated high doubles of the imaginary parts 
 *                of the backward products;
 *   backwardimlo are the updated low doubles of the imaginary parts 
 *                of the backward products;
 *   crossrehi    are the updated high doubles of the real parts
 *                of the cross products;
 *   crossrelo    are the updated low doubles of the real parts
 *                of the cross products;
 *   crossimhi    are the updated high doubles of the imaginary parts
 *                of the cross products;
 *   crossimlo    are the updated low doubles of the imaginary parts
 *                of the cross products. */

void CPU_dbl2_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   double ***forwardhi, double ***forwardlo,
   double ***backwardhi, double ***backwardlo,
   double ***crosshi, double ***crosslo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output,
 *   on real data.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csthi      high parts of constant coefficient series;
 *   cstlo      low parts of constant coefficient series;
 *   cffhi      high parts of coefficients of the monomials;
 *   cfflo      low parts of coefficients of the monomials;
 *   forwardhi  are all high parts of computed forward products,
 *   forwardlo  are all low parts of computed forward products,
 *   backwardhi are all high parts of computed backward products;
 *   backwardlo are all low parts of computed backward products;
 *   crosshi    are all high parts of computed cross products;
 *   crosslo    are all low parts of computed cross products.
 *
 * ON RETURN :
 *   outputhi   high parts of the values and all derivatives;
 *   outputlo   low parts of the values and all derivatives. */

void CPU_cmplx2_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   double ***forwardrehi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardrelo,
   double ***backwardimhi, double ***backwardimlo,
   double ***crossrehi, double ***crossrelo,
   double ***crossimhi, double ***crossimlo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cstrehi      high doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low doubles of the real parts
 *                of the constant coefficient series;
 *   cstimhi      high doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   forwardrehi  are all computed high doubles of the real parts
 *                of forward products;
 *   forwardrelo  are all computed low doubles of the real parts
 *                of forward products;
 *   forwardimhi  are all computed high doubles of the imaginary parts
 *                of forward products;
 *   forwardimlo  are all computed low doubles of the imaginary parts
 *                of forward products;
 *   backwardrehi are all computed high doubles of the real parts
 *                of the backward products;
 *   backwardrelo are all computed low doubles of the real parts
 *                of the backward products;
 *   backwardimhi are all computed high doubles of the imaginary parts
 *                of the backward products;
 *   backwardimlo are all computed low doubles of the imaginary parts
 *                of the backward products;
 *   crossrehi    are all computed high doubles of the real parts
 *                of the cross products;
 *   crossrelo    are all computed low doubles of the real parts
 *                of the cross products;
 *   crossimhi    are all computed high doubles of the imaginary parts
 *                of the cross products;
 *   crossimlo    are all computed low doubles of the imaginary parts
 *                of the cross products.
 *
 * ON RETURN :
 *   outputrehi   are the high doubles of the real parts
 *                of the values and all derivatives;
 *   outputrelo   are the low doubles of the real parts
 *                of the values and all derivatives;
 *   outputimhi   are the high doubles of the imaginary parts
 *                of the values and all derivatives;
 *   outputimlo   are the low doubles of the imaginary parts
 *                of the values and all derivatives. */

void CPU_dbl2_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   double ***forwardhi, double ***forwardlo,
   double ***backwardhi, double ***backwardlo,
   double ***crosshi, double ***crosslo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs, on real data.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csthi      high parts of constant coefficient series;
 *   cstlo      low parts of constant coefficient series;
 *   cffhi      high parts of coefficients of the monomials;
 *   cfflo      low parts of coefficients of the monomials;
 *   forwardhi  are all high parts of computed forward products,
 *   forwardlo  are all low parts of computed forward products,
 *   backwardhi are all high parts of computed backward products;
 *   backwardlo are all low parts of computed backward products;
 *   crosshi    are all high parts of computed cross products;
 *   crosslo    are all low parts of computed cross products.
 *   jobs       defines the addition jobs;
 *   verbose    if true, then output is written.
 *
 * ON RETURN :
 *   outputhi   high parts of the values and all derivatives;
 *   outputlo   low parts of the values and all derivatives. */

void CPU_cmplx2_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   double ***forwardrehi, double ***forwardrelo,
   double ***forwardimhi, double ***forwardimlo,
   double ***backwardrehi, double ***backwardrelo,
   double ***backwardimhi, double ***backwardimlo,
   double ***crossrehi, double ***crossrelo,
   double ***crossimhi, double ***crossimlo,
   AdditionJobs jobs, bool verbose );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs, on complex data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cstrehi      high doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low doubles of the real parts
 *                of the constant coefficient series;
 *   cstimhi      high doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   forwardrehi  are all computed high doubles of the real parts
 *                of forward products;
 *   forwardrelo  are all computed low doubles of the real parts
 *                of forward products;
 *   forwardimhi  are all computed high doubles of the imaginary parts
 *                of forward products;
 *   forwardimlo  are all computed low doubles of the imaginary parts
 *                of forward products;
 *   backwardrehi are all computed high doubles of the real parts
 *                of the backward products;
 *   backwardrelo are all computed low doubles of the real parts
 *                of the backward products;
 *   backwardimhi are all computed high doubles of the imaginary parts
 *                of the backward products;
 *   backwardimlo are all computed low doubles of the imaginary parts
 *                of the backward products;
 *   crossrehi    are all computed high doubles of the real parts
 *                of the cross products;
 *   crossrelo    are all computed low doubles of the real parts
 *                of the cross products;
 *   crossimhi    are all computed high doubles of the imaginary parts
 *                of the cross products;
 *   crossimlo    are all computed low doubles of the imaginary parts
 *                of the cross products.
 *   jobs         defines the addition jobs;
 *   verbose      if true, then output is written.
 *
 * ON RETURN :
 *   outputrehi   are the high doubles of the real parts
 *                of the values and all derivatives;
 *   outputrelo   are the low doubles of the real parts
 *                of the values and all derivatives;
 *   outputimhi   are the high doubles of the imaginary parts
 *                of the values and all derivatives;
 *   outputimlo   are the low doubles of the imaginary parts
 *                of the values and all derivatives. */

void CPU_dbl2_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo, double **outputhi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   on real data.
 *
 * ON ENTRY :
 *   dim        total number of variables;
 *   nbr        number of monomials, excluding the constant term;
 *   deg        truncation degree of the series;
 *   nvr        nvr[k] holds the number of variables in monomial k;
 *   idx        idx[k] has as many indices as the value of nvr[k],
 *              idx[k][i] defines the place of the i-th variable,
 *              with input values in input[idx[k][i]];
 *   csthi      high parts of constant coefficient series;
 *   cstlo      low parts of constant coefficient series;
 *   cffhi      high parts of coefficients of the monomials;
 *   cfflo      low parts of coefficients of the monomials;
 *   inputhi    has the high parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   inputlo    has the low parts of coefficients of the power series
 *              for all variables in the polynomial;
 *   outputhi   has space for high parts of the value and all derivatives;
 *   outputlo   has space for low parts of the value and all derivatives;
 *   cnvjobs    convolution jobs organized in layers;
 *   addjobs    addition jobs organized in layers;
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   outputhi   has the high parts of derivatives and the value,
 *              outputhi[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputhi[dim] contains the value of the polynomial;
 *   outputlo   has the low parts of derivatives and the value,
 *              outputlo[k], for k from 0 to dim-1, contains the
 *              derivative with respect to the variable k;
 *              outputlo[dim] contains the value of the polynomial;
 *   elapsedsec is the elapsed time in seconds. */

void CPU_cmplx2_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cstrehi      high doubles of the real parts
 *                of the constant coefficient series;
 *   cstrelo      low doubles of the real parts
 *                of the constant coefficient series;
 *   cstimhi      high doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cstimlo      low doubles of the imaginary parts
 *                of the constant coefficient series;
 *   cffrehi      cffrehi[k] has the deg+1 high doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffrelo      cffrelo[k] has the deg+1 low doubles of the real
 *                parts of the coefficient series of monomial k;
 *   cffimhi      cffimhi[k] has the deg+1 high doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   cffimlo      cffimlo[k] has the deg+1 low doubles of the imaginary
 *                parts of the coefficient series of monomial k;
 *   inputrehi    has the high doubles of the real parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputrelo    has the low doubles of the real part
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimhi    has the high doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   inputimlo    has the low doubles of the imaginary parts
 *                of the coefficients of the power series
 *                for all variables in the polynomial;
 *   outputrehi   has space for the high doubles of the real parts
 *                of the value and all derivatives;
 *   outputrelo   has space for the low doubles of the real parts
 *                of the value and all derivatives;
 *   outputimhi   has space for the high doubles of the imaginary parts
 *                of the value and all derivatives;
 *   outputimlo   has space for the low doubles of the imaginary parts
 *                of the value and all derivatives;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   vrblvl       is the verbose level.
 *
 * ON RETURN :
 *   outputrehi   has the high doubles of the real parts
 *                of the value and the derivatives,
 *   outputrelo   has the low doubles of the real parts
 *                of the value and the derivatives,
 *   outputimhi   has the high doubles of the imaginary parts
 *                of the value and the derivatives,
 *   outputimlo   has the low doubles of the imaginary parts
 *                of the value and the derivatives,
 *                output[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                output[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

#endif
