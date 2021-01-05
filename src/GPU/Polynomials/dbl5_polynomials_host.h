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
