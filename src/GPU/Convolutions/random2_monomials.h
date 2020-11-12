// The file random2_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in double double precision.

#ifndef __random2_monomials_h__
#define __random2_monomials_h__

bool make_real2_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhi, double *cfflo );
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
 *   cffhi   space allocated for deg+1 doubles;
 *   cfflo   space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffhi   deg+1 doubles with the high doubles of the coefficients;
 *   cfflo   deg+1 doubles with the low doubles of the coefficients. */

bool make_complex2_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo );
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
 *   cffrehi has space allocated for deg+1 doubles,
 *           for the high doubles of the real parts of the coefficients;
 *   cffrelo has space allocated for deg+1 doubles,
 *           for the low doubles of the real parts of the coefficients;
 *   cffimhi has space allocated for deg+1 doubles,
 *           for the high doubles of the imaginary parts of the coefficients;
 *   cffimlo has space allocated for deg+1 doubles,
 *           for the low doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffrehi holds deg+1 doubles with the high doubles of the real parts
 *           of the coefficients of the power series;
 *   cffrelo holds deg+1 doubles with the low doubles of the real parts
 *           of the coefficients of the power series;
 *   cffimhi holds deg+1 doubles with the high doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimlo holds deg+1 doubles with the low doubles of the imaginary parts
 *           of the coefficients of the power series. */

void make_real2_input ( int dim, int deg, double **datahi, double **datalo );
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
 *   datahi   space allocated for dim arrays of deg+1 doubles;
 *   datalo   space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datahi   datahi[i][j] is the high double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalo   datalo[i][j] is the low double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg. */

void make_complex2_input
 ( int dim, int deg, double **datarehi, double **datarelo,
   double **dataimhi, double **dataimlo );
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
 *   datarehi has space allocated for the high doubles of the real parts
 *            of dim series of degree deg.
 *   datarelo has space allocated for the low doubles of the real parts
 *            of dim series of degree deg;
 *   dataimhi has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg;
 *   dataimlo has space allocated for the low doubles of the imaginary parts
 *            of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehi contains the high doubles of the real parts of the data,
 *   datarehi contains the low doubles of the real parts of the data,
 *   dataimhi contains the high doubles of the imaginary parts of the data,
 *   dataimlo contains the low doubles of the imaginary parts of the data,
 *            data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

#endif
