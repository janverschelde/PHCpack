/* The file dbl4_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in quad double precision.
 *
 * The algorithmic differentiation is organized in two ways:
 * (1) CPU_dbl4_poly_evaldiff serves to verify the correctness;
 * (2) CPU_dbl4_poly_evaldiffjobs prepares the accelerated version,
 * with layered convolution jobs. */

#ifndef __dbl4_polynomials_host_h__
#define __dbl4_polynomials_host_h__

#include "convolution_jobs.h"
#include "addition_jobs.h"

void CPU_dbl4_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo, bool verbose=false );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in quad double precision.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihi    has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohi    has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilo    has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolo    has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
 *   forwardhihi  is work space for the highest doubles of nvr forward 
 *                products, forwardhihi[k] can store deg+1 doubles;
 *   forwardlohi  is work space for the second highest doubles of nvr
 *                forward  products, forwardlohi[k] can store deg+1 doubles;
 *   forwardhilo  is work space for the second lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   forwardlolo  is work space for the lowest doubles of nvr
 *                forward products, forwardlolo[k] can store deg+1 doubles;
 *   backwardhihi is work space for the highest doubles of nvr-1 backward
 *                products, backwardhihi[k] can store deg+1 doubles;
 *   backwardlohi is work space for the second highest doubles of nvr-1
 *                backward products, backwardlohi[k] can store deg+1 doubles;
 *   backwardhilo is work space for the second lowest doubles of nvr-1
 *                backward products, backwardhilo[k] can store deg+1 doubles;
 *   backwardlolo is work space for the lowest doubles of nvr- backward
 *                products, backwardlolo[k] can store deg+1 doubles;
 *   crosshihi    is work space for the highest doubles of nvr-2 cross
 *                products, crosshihi[k] can store deg+1 doubles;
 *   crosslohi    is work space for the second highest doubles of nvr-2
 *                cross products, crosslohi[k] can store deg+1 doubles;
 *   crosshilo    is work space for the second lowest doubles of nvr-2
 *                cross products, crosshilo[k] can store for deg+1 doubles.
 *   crosslolo    is work space for the lowest doubles of nvr-2 cross
 *                products, crosslolo[k] can store for deg+1 doubles.
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial. */

void CPU_dbl4_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihi    has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohi    has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilo    has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolo    has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

void CPU_dbl4_conv_job
 ( int deg, int nvr, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **forwardhihi, double **forwardlohi,
   double **forwardhilo, double **forwardlolo,
   double **backwardhihi, double **backwardlohi,
   double **backwardhilo, double **backwardlolo,
   double **crosshihi, double **crosslohi,
   double **crosshilo, double **crosslolo,
   ConvolutionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Computes one convolution defined by the job.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   nvr          number of variables in the monomial;
 *   idx          indices to the variables in the monomial;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihi    has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohi    has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilo    has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolo    has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
 *   forwardhihi  is work space for the highest doubles of nvr forward 
 *                products, forwardhihi[k] can store deg+1 doubles;
 *   forwardlohi  is work space for the second highest doubles of nvr
 *                forward  products, forwardlohi[k] can store deg+1 doubles;
 *   forwardhilo  is work space for the second lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   forwardlolo  is work space for the lowest doubles of nvr
 *                forward products, forwardhilo[k] can store deg+1 doubles;
 *   backwardhihi is work space for the highest doubles of nvr-1 backward
 *                products, backwardhihi[k] can store deg+1 doubles;
 *   backwardlohi is work space for the second highest doubles of nvr-1
 *                backward products, backwardlohi[k] can store deg+1 doubles;
 *   backwardhilo is work space for the second lowest doubles of nvr-1
 *                backward products, backwardhilo[k] can store deg+1 doubles;
 *   backwardlolo is work space for the lowest doubles of nvr-1 backward
 *                products, backwardlolo[k] can store deg+1 doubles;
 *   crosshihi    is work space for the highest doubles of nvr-2 cross
 *                products, crosshihi[k] can store deg+1 doubles;
 *   crosslohi    is work space for the second highest doubles of nvr-2
 *                cross products, crosslohi[k] can store deg+1 doubles;
 *   crosshilo    is work space for the second lowest doubles of nvr-2
 *                cross products, crosshilo[k] can store for deg+1 doubles;
 *   crosslolo    is work space for the lowest doubles of nvr-2 cross
 *                products, crosslolo[k] can store for deg+1 doubles;
 *   job          defines the convolution job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihi  are the updated highest parts of forward products;
 *   forwardlohi  are the updated second highest parts of forward products;
 *   forwardhilo  are the updated second lowest parts of forward products;
 *   forwardlolo  are the updated lowest parts forward products;
 *   backwardhihi are the updated highest parts of backward products;
 *   backwardlohi are the updated second highest parts of backward products;
 *   backwardhilo are the updated second lowest parts of backward products;
 *   backwardlolo are the updated lowest parts backward products;
 *   crosshihi    are the updated highest parts of cross products;
 *   crosslohi    are the updated second highest parts of cross products;
 *   crosshilo    are the updated second lowest parts of cross products;
 *   crosslolo    are the updated lowest parts cross products. */

void CPU_dbl4_add_job
 ( int deg,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJob job, bool verbose );
/*
 * DESCRIPTION :
 *   Does one update defined by the job.
 *
 * ON ENTRY :
 *   deg          degree of the series;
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   forwardhihi  are all highest parts of computed forward products,
 *   forwardlohi  are all second highest parts of computed forward products,
 *   forwardhilo  are all second lowest parts of computed forward products,
 *   forwardlolo  are all lowest parts of computed forward products,
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products;
 *   job          defines the addition job;
 *   verbose      if true, then is verbose.
 *
 * ON RETURN :
 *   forwardhihi  are the updated highest parts of forward products;
 *   forwardlohi  are the updated second highest parts of forward products;
 *   forwardhilo  are the updated second lowest parts of forward products;
 *   forwardlolo  are the updated lowest parts forward products;
 *   backwardhihi are the updated highest parts of backward products;
 *   backwardlohi are the updated second highest parts of backward products;
 *   backwardhilo are the updated second lowest parts of backward products;
 *   backwardlolo are the updated lowest parts backward products;
 *   crosshihi    are the updated highest parts of cross products;
 *   crosslohi    are the updated second highest parts of cross products;
 *   crosshilo    are the updated second lowest parts of cross products;
 *   crosslolo    are the updated lowest parts cross products. */

void CPU_dbl4_poly_updates
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions in a straightforward manner to the final output.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   forwardhihi  are all highest parts of computed forward products;
 *   forwardlohi  are all second highest parts of computed forward products;
 *   forwardhilo  are all second lowest parts of computed forward products;
 *   forwardlolo  are all lowest parts of computed forward products;
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products.
 *
 * ON RETURN :
 *   outputhihi   highest parts of the values and all derivatives;
 *   outputlohi   second highest parts of the values and all derivatives;
 *   outputhilo   second lowest parts of the values and all derivatives;
 *   outputlolo   lowest parts of the values and all derivatives. */

void CPU_dbl4_poly_addjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   double ***forwardhihi, double ***forwardlohi,
   double ***forwardhilo, double ***forwardlolo,
   double ***backwardhihi, double ***backwardlohi,
   double ***backwardhilo, double ***backwardlolo, 
   double ***crosshihi, double ***crosslohi,
   double ***crosshilo, double ***crosslolo,
   AdditionJobs jobs, bool verbose=false );
/*
 * DESCRIPTION :
 *   Given the forward, backward, and cross products for every monomial,
 *   makes all additions as defined by the addition jobs.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   forwardhihi  are all highest parts of computed forward products;
 *   forwardlohi  are all second highest parts of computed forward products;
 *   forwardhilo  are all second lowest parts of computed forward products;
 *   forwardlolo  are all lowest parts of computed forward products;
 *   backwardhihi are all highest parts of computed backward products;
 *   backwardlohi are all second highest parts of computed backward products;
 *   backwardhilo are all second lowest parts of computed backward products;
 *   backwardlolo are all lowest parts of computed backward products;
 *   crosshihi    are all highest parts of computed cross products;
 *   crosslohi    are all second highest parts of computed cross products;
 *   crosshilo    are all second lowest parts of computed cross products;
 *   crosslolo    are all lowest parts of computed cross products;
 *   jobs         defines the addition jobs;
 *   verbose      if true, then output is written.
 *
 * ON RETURN :
 *   outputhihi   highest parts of the values and all derivatives;
 *   outputlohi   second highest parts of the values and all derivatives;
 *   outputhilo   second lowest parts of the values and all derivatives;
 *   outputlolo   lowest parts of the values and all derivatives. */

void CPU_dbl4_poly_evaldiffjobs
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo, 
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *elapsedsec, bool verbose=false );
/*
 * DESCRIPTION :
 *   Computes the convolutions in the order as defined by cnvjobs,
 *   performs the updates to the values as defined by addjobs,
 *   all other parameters are the same as in the other function.
 *
 * ON ENTRY :
 *   dim          total number of variables;
 *   nbr          number of monomials, excluding the constant term;
 *   deg          truncation degree of the series;
 *   nvr          nvr[k] holds the number of variables in monomial k;
 *   idx          idx[k] has as many indices as the value of nvr[k],
 *                idx[k][i] defines the place of the i-th variable,
 *                with input values in input[idx[k][i]];
 *   csthihi      highest parts of constant coefficient series;
 *   cstlohi      second highest parts of constant coefficient series;
 *   csthilo      second lowest parts of constant coefficient series;
 *   cstlolo      lowest parts of constant coefficient series;
 *   cffhihi      cffhihi[k] has deg+1 doubles for the highest parts
 *                of the coefficient series of monomial k;
 *   cfflohi      cfflohi[k] has deg+1 doubles for the second highest parts
 *                of the coefficient series of monomial k;
 *   cffhilo      cffhilo[k] has deg+1 doubles for the second lowest parts
 *                of the coefficient series of monomial k;
 *   cfflolo      cfflolo[k] has deg+1 doubles for the lowest parts
 *                of the coefficient series of monomial k;
 *   inputhihi    has the highest parts of the power series
 *                for all variables in the polynomial;
 *   inputlohi    has the second highest parts of the power series
 *                for all variables in the polynomial;
 *   inputhilo    has the second lowest parts of the power series
 *                for all variables in the polynomial;
 *   inputlolo    has the lowest parts of the power series
 *                for all variables in the polynomial;
 *   outputhihi   has space allocated for dim+1 series of degree deg;
 *   outputlohi   has space allocated for dim+1 series of degree deg;
 *   outputhilo   has space allocated for dim+1 series of degree deg;
 *   outputlolo   has space allocated for dim+1 series of degree deg;
 *   cnvjobs      convolution jobs organized in layers;
 *   addjobs      addition jobs organized in layers;
 *   verbose      if true, writes one line to screen for every convolution.
 *
 * ON RETURN :
 *   outputhihi   has the highest parts of derivatives and the value,
 *                outputhihi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhihi[dim] contains the value of the polynomial;
 *   outputlohi   has the second highest parts of derivatives and the value,
 *                outputlohi[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlohi[dim] contains the value of the polynomial;
 *   outputhilo   has the second lowest parts of derivatives and the value,
 *                outputhilo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputhilo[dim] contains the value of the polynomial;
 *   outputlolo   has the lowest parts of derivatives and the value,
 *                outputlolo[k], for k from 0 to dim-1, contains the
 *                derivative with respect to the variable k;
 *                outputlolo[dim] contains the value of the polynomial;
 *   elapsedsec   is the elapsed time in seconds. */

#endif
