// The file dbl_systems_kernels.h specifies functions to define
// memory transfers ad kernel launches to evaluate and differentiate
// monomials with common factors in double precision.

#ifndef __dbl_systems_kernels_h__
#define __dbl_systems_kernels_h__

#include "convolution_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"

void write_arithmetic_intensity
 ( long long int flopcnt, long long int bytecnt,
   double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the arithmetic intensity, the kernel and wall flops.
 *
 * ON ENTRY :
 *   flopcnt  counts of the flops;
 *   bytecnt  counts of the bytes;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void write_dbl_cnvflops
 ( int dim, int deg, int ctype,
   ConvolutionJobs cnvjobs, double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the kernel and wall time flops for the convolution jobs.
 *
 * ON ENTRY :
 *   dim      number of series;
 *   deg      order of the series, truncation degree;
 *   ctype    0 if on real data, 1 if on complex data;
 *   cnvjobs  defines the convolution jobs;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void write_vectorized_cnvincflops
 ( int dim, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double kernms, double wallsec );
/*
 * DESCRIPTION :
 *   Writes the kernel and wall time flops for the convolution
 *   and the increment jobs in the complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim      number of series;
 *   deg      order of the series, truncation degree;
 *   cnvjobs  defines the convolution jobs;
 *   incjobs  defines the increment jobs;
 *   kernms   kernel time elapsed in milliseconds; 
 *   wallsec  wall clock time elapsed in seconds. */

void dbl_evaldiffdata_to_output
 ( double *data, double ***output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the real data computed on the device to the output.
 *
 * ON ENTRY :
 *   data     coefficients of all monomials and input series, 
 *            computed forward, backward, and cross products;
 *   output   space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   output   output[i][dim] contains the power series value
 *            of the i-th monomial, and
 *            output[i][k] contains the power series value of
 *            the k-th derivative of the i-th monomial. */

void cmplx_evaldiffdata_to_output
 ( double *datare, double *dataim, double ***outputre, double ***outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device to the output.
 *
 * ON ENTRY :
 *   datare   real parts of coefficients of all monomials and input, 
 *            computed forward, backward, and cross products;
 *   dataim   imaginary parts of coefficients of all monomials and input, 
 *            computed forward, backward, and cross products;
 *   outputre has space for the value and all derivatives;
 *   outputim has space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputre, in outputre[i][dim] are the real parts of the value
 *            of the i-th monomial, and
 *            outputre[i][k] has the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputim, in outputre[i][dim] are the imaginary parts of the value
 *            of the i-th monomial, and
 *            outputim[i][k] has the imaginary parts of the value of
 *            the k-th derivative of the i-th monomial. */

void cmplxvectorized_evaldiffdata_to_output
 ( double *datari, double ***outputre, double ***outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totalcff, int offsetri, int vrblvl );
/*
 * DESCRIPTION :
 *   Extracts the complex data computed on the device
 *   using vectorized arithmetic to the output.
 *
 * ON ENTRY :
 *   datari   coefficients of all monomials and input, 
 *            computed forward, backward, and cross products;
 *   outputre has space for the value and all derivatives;
 *   outputim has space for the value and all derivatives;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] is the number of variables for monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   fstart   fstart[k] has the start position of the forward products
 *            for the k-th monomial;
 *   bstart   fstart[k] has the start position of the backward products
 *            for the k-th monomial;
 *   cstart   fstart[k] has the start position of the cross products
 *            for the k-th monomial;
 *   totalcff is the total number of coefficients without vectorization;
 *   offsetri is the size of the second operand;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputre, in outputre[i][dim] are the real parts of the value
 *            of the i-th monomial, and
 *            outputre[i][k] has the real parts of the value of
 *            the k-th derivative of the i-th monomial;
 *   outputim, in outputre[i][dim] are the imaginary parts of the value
 *            of the i-th monomial, and
 *            outputim[i][k] has the imaginary parts of the value of
 *            the k-th derivative of the i-th monomial. */

void GPU_dbl_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cff, double **input, double ***output,
   ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs, on real data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cff      cff[k] has deg+1 doubles for the coefficient series
 *            of monomial k;
 *   input    contains the coefficients of the power series
 *            for all variables in the polynomial;
 *   output   space allocated for the value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplx_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffre, double **cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim, ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   on complex data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cffre    cffre[k] has deg+1 doubles of the real parts for the
 *            coefficient series of monomial k;
 *   cffim    cffim[k] has deg+1 doubles of the imaginary parts for the
 *            coefficient series of monomial k;
 *   inputre  real parts of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputim  imaginary parts of the coefficients of the series
 *            for all variables in the polynomial;
 *   outputre has space for real parts of value and all derivatives;
 *   outputim has space for imaginary parts of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputre has real parts of derivatives and the value,
 *            outputre[k], for k from 0 to dim-1, has the real part
 *            of the derivative with respect to the variable k;
 *            outputre[dim] contains the real part of the value;
 *   outputim has imaginary parts of derivatives and the value,
 *            outputim[k], for k from 0 to dim-1, has the imaginary part
 *            of the derivative with respect to the variable k;
 *            outputim[dim] contains the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_cmplxvectorized_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffre, double **cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates and differentiates a monomial system.
 *   Computes the convolutions in the order as defined by jobs,
 *   on complex data using complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   szt      number of threads in a block, must equal deg + 1;
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant term;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cffre    cffre[k] has deg+1 doubles of the real parts for the
 *            coefficient series of monomial k;
 *   cffim    cffim[k] has deg+1 doubles of the imaginary parts for the
 *            coefficient series of monomial k;
 *   inputre  real parts of the coefficients of the series
 *            for all variables in the polynomial;
 *   inputim  imaginary parts of the coefficients of the series
 *            for all variables in the polynomial;
 *   outputre has space for real parts of value and all derivatives;
 *   outputim has space for imaginary parts of value and all derivatives;
 *   cnvjobs  convolution jobs organized in layers;
 *   incjobs  increment jobs organized in layers;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   outputre has real parts of derivatives and the value,
 *            outputre[k], for k from 0 to dim-1, has the real part
 *            of the derivative with respect to the variable k;
 *            outputre[dim] contains the real part of the value;
 *   outputim has imaginary parts of derivatives and the value,
 *            outputim[k], for k from 0 to dim-1, has the imaginary part
 *            of the derivative with respect to the variable k;
 *            outputim[dim] contains the imaginary part of the value;
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

void GPU_dbl_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cff, double *acc, double **input, double ***output,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cff       coefficients of the monomials;
 *   acc       space to accumulate one power series of degree deg;
 *   input     coefficients of the power series of degree deg,
 *             for dim variables;
 *   output    space for the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl   is the verbose level, if zero, then no output.
 *
 * ON RETURN :
 *   cff       contains the evaluated common factors;
 *   output    evaluated and differentiated monomials in the system,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             output[i][idx[k]] is the derivative w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffre     real parts of coefficients of the monomials;
 *   cffim     imaginary parts of coefficients of the monomials;
 *   accre     space to accumulate one power series of degree deg;
 *   accim     space to accumulate one power series of degree deg;
 *   inputre   real parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   inputim   imaginary parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   outputre  has space for the real parts of the output;
 *   outputim  has space for the imaginary parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffre     has the real parts of the evaluated common factors;
 *   cffim     has the imaginary parts of the evaluated common factors;
 *   outputre  real parts of evaluated and differentiated monomials,
 *             outputre[i][dim] is the real part of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputre[i][idx[k]] is the real part of derivative
 *             w.r.t. idx[k];
 *   outputim  imaginary parts of evaluated and differentiated monomials,
 *             outputim[i][dim] is the imaginary part of the value of the
 *             input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputim[i][idx[k]] is the imaginary part of te derivative
 *             w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplxvectorized_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
   double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series using complex vectorized arithmetic.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffre     real parts of coefficients of the monomials;
 *   cffim     imaginary parts of coefficients of the monomials;
 *   accre     space to accumulate one power series of degree deg;
 *   accim     space to accumulate one power series of degree deg;
 *   inputre   real parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   inputim   imaginary parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   outputre  has space for the real parts of the output;
 *   outputim  has space for the imaginary parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffre     has the real parts of the evaluated common factors;
 *   cffim     has the imaginary parts of the evaluated common factors;
 *   outputre  real parts of evaluated and differentiated monomials,
 *             outputre[i][dim] is the real part of the value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputre[i][idx[k]] is the real part of derivative
 *             w.r.t. idx[k];
 *   outputim  imaginary parts of evaluated and differentiated monomials,
 *             outputim[i][dim] is the imaginary part of the value of the
 *             input at the i-th monomial, and for k in range 0..nvr[i]-1:
 *             outputim[i][idx[k]] is the imaginary part of te derivative
 *             w.r.t. idx[k];
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_dbl_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx,
   double ***cff, double **input, double ***output,
   double **funval, double ***jacval, double *totcnvlapsedms, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cff       cff[i][j] is the coefficient of the j-th monomial
 *             in the i-th column;
 *   input     coefficients of the power series of degree deg,
 *             for dim variables;
 *   output    space for the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   output    used as work space for one column;
 *   funval    the evaluated series for each polynomial;
 *   jacval    matrix series of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

void GPU_cmplx_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffre, double ***cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim, 
   double **funvalre, double **funvalim,
   double ***jacvalre, double ***jacvalim, double *totcnvlapsedms,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cffre     cffre[i][j] are the real parts of the coefficients
 *             of the j-th monomial in the i-th column;
 *   cffim     cffim[i][j] are the imaginary parts of the coefficients
 *             of the j-th monomial in the i-th column;
 *   inputre   real parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   inputim   imaginary parts of coefficients of the series
 *             of degree deg, for dim variables;
 *   outputre  has space for the real parts of the output;
 *   outputim  has space for the imaginary parts of the output;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   outputre  real parts of evaluated and differentiated monomials,
 *             used as work space for one column;
 *   outputim  imaginary parts of evaluated and differentiated monomials,
 *             used as work space for one column;
 *   funvalre  real parts of the evaluated series for each polynomial;
 *   funvalim  imaginary parts of the evaluated series for each polynomial;
 *   jacvalre  real parts of the matrix series of all derivatives;
 *   jacvalim  imaginary parts of the matrix series of all derivatives;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions. */

#endif
