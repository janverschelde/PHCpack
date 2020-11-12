// The file random4_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in quad double precision.

#ifndef __random4_monomials_h__
#define __random4_monomials_h__

bool make_real4_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo );
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
 *   cffhihi has space allocated for deg+1 doubles;
 *   cfflohi has space allocated for deg+1 doubles;
 *   cffhilo has space allocated for deg+1 doubles;
 *   cfflolo has space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffhihi stores the deg+1 highest doubles of the coefficients;
 *   cfflohi stores the deg+1 second highest doubles of the coefficients;
 *   cffhilo stores the deg+1 second lowest doubles of the coefficients;
 *   cfflolo stores the deg+1 lowest doubles of the coefficients. */

bool make_complex4_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi,
   double *cffimhilo, double *cffimlolo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nvr       number of variables with positive power in the monomial;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   exp       space allocated for nvr integers;
 *   cffrehihi has space allocated for deg+1 doubles,
 *             for the highest doubles of the real parts of the coefficients;
 *   cffrelohi has space allocated for deg+1 doubles, for the second
 *             highest doubles of the real parts of the coefficients;
 *   cffrehilo has space allocated for deg+1 doubles, for the second
 *             lowest doubles of the real parts of the coefficients;
 *   cffrelolo has space allocated for deg+1 doubles, for the lowest
 *             doubles of the real parts of the coefficients;
 *   cffimhihi has space allocated for deg+1 doubles, for the highest
 *             doubles of the imaginary parts of the coefficients;
 *   cffimlohi has space allocated for deg+1 doubles, for the second
 *             highest doubles of the imaginary parts of the coefficients;
 *   cffimhilo has space allocated for deg+1 doubles, for the second
 *             lowest doubles of the imaginary parts of the coefficients.
 *   cffimlolo has space allocated for deg+1 doubles, for the lowest
 *             doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx       nvr integers in the range from 0 to dim-1,
 *             idx(k) is the index of the k-th variable in the monomial;
 *   exp       nvr positive integers with the powers of the variables,
 *             exp(k) is the power of the variable with index idx(k);
 *   cffrehihi holds the deg+1 highest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrelohi holds the deg+1 second highest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrehilo holds the deg+1 second lowest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffrelolo holds the deg+1 lowest doubles of the real parts
 *             of the coefficients of the power series;
 *   cffimhihi holds the deg+1 highest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimlohi holds the deg+1 second highest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimhilo holds the deg+1 second lowest doubles of the imaginary parts
 *             of the coefficients of the power series;
 *   cffimlolo holds deg+1 lowest doubles of the imaginary parts
 *             of the coefficients of the power series. */

void make_real4_input
 ( int dim, int deg, double **datahihi, double **datalohi,
                     double **datahilo, double **datalolo );
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
 *   datahihi has space allocated for dim arrays of deg+1 doubles;
 *   datalohi has space allocated for dim arrays of deg+1 doubles;
 *   datahilo has space allocated for dim arrays of deg+1 doubles;
 *   datalolo has space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datahihi holds the highest doubles of the input series,
 *            datahihi[i][j] is the highest double of the j-th coefficient 
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalohi holds the second highest doubles of the input series,
 *            datalohi[i][j] is the 2nd highest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datahilo holds the second lowest doubles of the input series,
 *            datalolo[i][j] is the 2nd lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalolo holds the lowest doubles of the input series,
 *            datalolo[i][j] is the lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg. */

void make_complex4_input
 ( int dim, int deg,
   double **datarehihi, double **datarelohi,
   double **datarehilo, double **datarelolo,
   double **dataimhihi, double **dataimlohi,
   double **dataimhilo, double **dataimlolo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim        dimension of the input;
 *   deg        degree of the power series;
 *   datarehihi has space allocated for the highest doubles
 *              of the real parts of dim series of degree deg.
 *   datarelohi has space allocated for the second highest doubles
 *              of the real parts of dim series of degree deg.
 *   datarehilo has space allocated for the second lowest doubles
 *              of the real parts of dim series of degree deg;
 *   datarelolo has space allocated for the lowest doubles
 *              of the real parts of dim series of degree deg;
 *   dataimhihi has space allocated for the highest doubles
 *              of the imaginary parts of dim series of degree deg;
 *   dataimlohi has space allocated for the second highest doubles
 *              of the imaginary parts of dim series of degree deg;
 *   dataimhilo has space allocated for the second lowest doubles
 *              of the imaginary parts of dim series of degree deg.
 *   dataimlolo has space allocated for the lowest doubles
 *              of the imaginary parts of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehihi stores the highest doubles of the real parts,
 *   datarelohi stores the second highest doubles of the real parts,
 *   datarehilo stores the second lowest doubles of the real parts,
 *   datarelolo stores the lowest doubles of the real parts,
 *   dataimhihi stores the highest doubles of the imaginary parts,
 *   dataimlohi stores the second highest doubles of the imaginary parts,
 *   dataimhilo stores the second lowest doubles of the imaginary parts,
 *   dataimlolo stores the lowest doubles of the imaginary parts,
 *              data[i][j] is the j-th coefficient of the i-th series,
 *              for i in 0..dim-1 and j in 0..deg. */

#endif
