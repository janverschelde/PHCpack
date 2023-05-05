// The file dbl_fabry_host.h specifies the functions for the Fabry ratios
// of series in double precision on real numbers.

#ifndef __dbl_fabry_host_h__
#define __dbl_fabry_host_h__

int dbl_fabry_ratio
 ( int deg, double *input, double *ratio, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the ratio of the next to last over the last coefficient
 *   in the input series.
 *
 * ON ENTRY :
 *   deg       degree of the power series;
 *   input     coefficients of a series of degree deg;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   ratio     contains the ratio for the series of the input,
 *             in absolute value. */

int dbl_fabry_vector
 ( int dim, int deg, double **input, double *ratios, int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the ratios of the next to last over the last coefficient
 *   in the input series, for all components.
 *
 * ON ENTRY :
 *   dim       number of series;
 *   deg       degree of the power series;
 *   input     coefficients of series of degree deg, for dim variables;
 *   ratios    has space for dim doubles;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   ratios    contains the ratio for each series of the input. */

int dbl_fabry_smallest
 ( int dim, double *ratios, double *step, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the Fabry ratios for every series,
 *   the Fabry step is the smallest of the ratios.
 *
 * ON ENTRY :
 *   dim       number of series, size of ratios;
 *   ratios    ratio for each series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   step      the smallest number in ratios. */

int dbl_fabry_step
 ( int dim, int deg, double **input, double *ratios, double *step,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Computes the Fabry ratios from the input series
 *   and the smallest ratio as the step.
 *
 * ON ENTRY :
 *   dim       number of series;
 *   deg       degree of the power series;
 *   input     coefficients of series of degree deg, for dim variables;
 *   ratios    has space for dim doubles;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   ratios    contains the ratio for each series of the input;
 *   step      is the smallest number in ratios. */

int dbl_fabry_predict1
 ( int deg, double *input, double step, double *output, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates a rational approximation with linear denominator
 *   at the step size to predict one coordinate of the next point.
 *   If a coefficient of the input series is too small for a division,
 *   then the series is evaluated and 1 is returned.
 *
 * ON ENTRY :
 *   deg       degree of the power series;
 *   input     coefficients of a series of degree deg;
 *   step      the step size is the value for the continuation parameter;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   output    predicted value of one coordinate of the next point. */

int dbl_fabry_predictor
 ( int dim, int deg, double **input, double step, double *output,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Given on input a sequence of series, predicts the next point,
 *   using a rational approximation with linear denominator.
 *   The integer on return is the number of failures, that is:
 *   the number of times a coefficient was too small for a division.
 *
 * ON ENTRY :
 *   dim       number of series;
 *   deg       degree of the power series;
 *   input     coefficients of series of degree deg, for dim variables;
 *   step      value for the continuation parameter;
 *   output    has space for dim doubles;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   output    prediction for the next point on a path. */

#endif
