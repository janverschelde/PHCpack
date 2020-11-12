// The file random_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in double precision.

#ifndef __random_monomials_h__
#define __random_monomials_h__

bool sorted_insert ( int n, int *data );
/*
 * DESCRIPTION :
 *   Inserts data[n] in the sequence of n sorted numbers in data.
 *   Returns true if data[n] was already inserted, that is:
 *   there is an index k less than n, for which data[k] == data[n].
 *   Returns false if data[n] is not a duplicate number. */

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp, double *cff );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cff     space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cff     deg+1 doubles with the coefficients of the power series. */

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffre, double *cffim );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cffre   space allocated for deg+1 doubles,
 *           for the real parts of the coefficients;
 *   cffim   space allocated for deg+1 doubles,
 *           for the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffre   deg+1 doubles with the real parts
 *           of the coefficients of the power series;
 *   cffim   deg+1 doubles with the imaginary parts
 *           of the coefficients of the power series. */

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac );
/*
 * DESCRIPTION :
 *   Extracts all exponents strictly larger than one.
 *  
 * ON ENTRY :
 *   nvr     number of variables in exp, exp[k] >= 1,
 *           for all k from 0 to nvr-1;
 *   exp     exponents of a monomial;
 *   expfac  space for nvr integers.
 *
 * ON RETURN :
 *   nbrfac  number of exponents in exp strictly larger than one;
 *   expfac  exponents of the common factor,
 *           if exp[k] > 1, then expfac[k] = exp[k]-1.  */

void make_real_input ( int dim, int deg, double **data );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   data     space allocated for dim series of degree deg.
 *
 * ON RETURN :
 *   data     data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

void make_complex_input
 ( int dim, int deg, double **datare, double **dataim );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datare   space allocated for the real parts of dim series
 *            of degree deg.
 *   dataim   space allocated for the imaginary parts of dim series
 *            of degree deg.
 *
 * ON RETURN :
 *   datare   the real parts of the data,
 *   dataim   the imaginary parts of the data,
 *            data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

#endif
