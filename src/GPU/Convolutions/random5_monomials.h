// The file random5_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in penta double precision.

#ifndef __random5_monomials_h__
#define __random5_monomials_h__

bool make_real5_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk );
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
 *   cfftb   space allocated for deg+1 doubles;
 *   cffix   space allocated for deg+1 doubles;
 *   cffmi   space allocated for deg+1 doubles;
 *   cffrg   space allocated for deg+1 doubles;
 *   cffpk   space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cfftb   deg+1 doubles with the highest coefficient doubles;
 *   cffix   deg+1 doubles with the second highest coefficient doubles;
 *   cffmi   deg+1 doubles with the middle coefficient doubles;
 *   cffrg   deg+1 doubles with the second lowest coefficient doubles;
 *   cffpk   deg+1 doubles with the lowest coefficient doubles. */

bool make_complex5_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk );
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
 *   cffretb has space allocated for deg+1 highest doubles,
 *           for the real parts of the coefficients;
 *   cffreix has space allocated for deg+1 second highest doubles,
 *           for the real parts of the coefficients;
 *   cffremi has space allocated for deg+1 middle doubles,
 *           for the real parts of the coefficients;
 *   cffrerg has space allocated for deg+1 second lowest doubles,
 *           for the real parts of the coefficients;
 *   cffrepk has space allocated for deg+1 lowest doubles,
 *           for the real parts of the coefficients;
 *   cffimtb has space allocated for deg+1 highest doubles,
 *           for the imaginary parts of the coefficients;
 *   cffimix has space allocated for deg+1 second highest doubles,
 *           for the imaginary parts of the coefficients;
 *   cffimmi has space allocated for deg+1 middle doubles,
 *           for the imaginary parts of the coefficients;
 *   cffimrg has space allocated for deg+1 lowest doubles,
 *           for the imaginary parts of the coefficients;
 *   cffimpk has space allocated for deg+1 second lowest doubles,
 *           for the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffretb holds deg+1 highest doubles of the real parts
 *           of the coefficients of the power series;
 *   cffreix holds deg+1 second highest doubles of the real parts
 *           of the coefficients of the power series;
 *   cffremi holds deg+1 middle doubles of the real parts
 *           of the coefficients of the power series;
 *   cffrerg holds deg+1 second lowest doubles of the real parts
 *           of the coefficients of the power series;
 *   cffrepk holds deg+1 lowest doubles of the real parts
 *           of the coefficients of the power series;
 *   cffimtb holds deg+1 highest doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimix holds deg+1 second highest doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimmi holds deg+1 middle doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimrg holds deg+1 second lowest doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimpk holds deg+1 lowest doubles of the imaginary parts
 *           of the coefficients of the power series. */

void make_real5_input
 ( int dim, int deg, double **datatb, double **dataix,
   double **datami, double **datarg, double **datapk );
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
 *   datatb   space allocated for dim arrays of deg+1 doubles;
 *   dataix   space allocated for dim arrays of deg+1 doubles;
 *   datami   space allocated for dim arrays of deg+1 doubles;
 *   datarg   space allocated for dim arrays of deg+1 doubles;
 *   datapk   space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datatb   datatb[i][j] is the highest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   dataix   dataix[i][j] is the second highest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datami   datami[i][j] is the middle double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datarg   datarg[i][j] is the second lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datapk   datapk[i][j] is the lowest double of the j-th coefficient
 *            of the i-th series, for i in 0..dim-1 and j in 0..deg. */

void make_complex5_input
 ( int dim, int deg,
   double **dataretb, double **datareix, double **dataremi,
   double **datarerg, double **datarepk,
   double **dataimtb, double **dataimix, double **dataimmi,
   double **dataimrg, double **dataimpk );
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
 *   dataretb has space allocated for the highest doubles of the real parts
 *            of dim series of degree deg.
 *   datareix has space allocated for the second highest doubles of the
 *            real parts of dim series of degree deg.
 *   dataremi has space allocated for the middle doubles of the real parts
 *            of dim series of degree deg.
 *   datarerg has space allocated for the second lowest doubles of the
 *            real parts of dim series of degree deg;
 *   datarepk has space allocated for the lowest doubles of the real parts
 *            of dim series of degree deg;
 *   dataimtb has space allocated for the highest doubles of the
 *            imaginary parts of dim series of degree deg;
 *   dataimix has space allocated for the second highest doubles of the
 *            imaginary parts of dim series of degree deg;
 *   dataimmi has space allocated for the middle doubles of the
 *            imaginary parts of dim series of degree deg;
 *   dataimrg has space allocated for the second lowest doubles of the
 *            imaginary parts of dim series of degree deg;
 *   dataimpk has space allocated for the lowest doubles of the
 *            imaginary parts of dim series of degree deg.
 *
 * ON RETURN :
 *   dataretb stores the highest doubles of the real parts of the data,
 *   datareix stores the second highest doubles of the real parts,
 *   dataremi stores the middle doubles of the real parts of the data,
 *   datarerg stores the second lowest doubles of the real parts,
 *   datarepk stores the lowest doubles of the real parts of the data,
 *   dataimtb stores the highest doubles of the imaginary parts of the data,
 *   dataimix stores the second highest doubles of the imaginary parts,
 *   dataimmi stores the middle doubles of the imaginary parts of the data,
 *   dataimrg stores the second lowest doubles of the imaginary parts,
 *   dataimpk stores the lowest doubles of the imaginary parts of the data,
 *            data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

#endif
