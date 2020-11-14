// The file random10_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in deca double precision.

#ifndef __random10_monomials_h__
#define __random10_monomials_h__

bool make_real10_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrtb, double *cffrix, double *cffrmi,
   double *cffrrg, double *cffrpk,
   double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables with positive power in the monomial;
 *   pwr      largest power of a variable;
 *   deg      degree of the power series coefficient;
 *   exp      space allocated for nvr integers;
 *   cffrtb   space allocated for deg+1 doubles;
 *   cffrix   space allocated for deg+1 doubles;
 *   cffrmi   space allocated for deg+1 doubles;
 *   cffrrg   space allocated for deg+1 doubles;
 *   cffrpk   space allocated for deg+1 doubles;
 *   cffltb   space allocated for deg+1 doubles;
 *   cfflix   space allocated for deg+1 doubles;
 *   cfflmi   space allocated for deg+1 doubles;
 *   cfflrg   space allocated for deg+1 doubles;
 *   cfflpk   space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx      nvr integers in the range from 0 to dim-1,
 *            idx(k) is the index of the k-th variable in the monomial;
 *   exp      nvr positive integers with the powers of the variables,
 *            exp(k) is the power of the variable with index idx(k);
 *   cffrtb   deg+1 doubles with the highest coefficient doubles;
 *   cffrix   deg+1 doubles with the second highest coefficient doubles;
 *   cffrmi   deg+1 doubles with the third highest coefficient doubles;
 *   cffrrg   deg+1 doubles with the fourth highest coefficient doubles;
 *   cffrpk   deg+1 doubles with the fifth highest coefficient doubles;
 *   cffltb   deg+1 doubles with the fifth lowest coefficient doubles;
 *   cfflix   deg+1 doubles with the fourth lowest coefficient doubles;
 *   cfflmi   deg+1 doubles with the third lowest coefficient doubles;
 *   cfflrg   deg+1 doubles with the second lowest coefficient doubles;
 *   cfflpk   deg+1 doubles with the lowest coefficient doubles. */

bool make_complex10_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrertb, double *cffrerix, double *cffrermi,
   double *cffrerrg, double *cffrerpk,
   double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi,
   double *cffimrrg, double *cffimrpk,
   double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk );
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
 *   cffrertb has space allocated for deg+1 highest doubles,
 *            for the real parts of the coefficients;
 *   cffrerix has space allocated for deg+1 second highest doubles,
 *            for the real parts of the coefficients;
 *   cffrermi has space allocated for deg+1 third highest doubles,
 *            for the real parts of the coefficients;
 *   cffrerrg has space allocated for deg+1 fourth highest doubles,
 *            for the real parts of the coefficients;
 *   cffrerpk has space allocated for deg+1 fifth highest doubles,
 *            for the real parts of the coefficients;
 *   cffreltb has space allocated for deg+1 fifth lowest doubles,
 *            for the real parts of the coefficients;
 *   cffrelix has space allocated for deg+1 fourth lowest doubles,
 *            for the real parts of the coefficients;
 *   cffrelmi has space allocated for deg+1 third lowest doubles,
 *            for the real parts of the coefficients;
 *   cffrelrg has space allocated for deg+1 second lowest doubles,
 *            for the real parts of the coefficients;
 *   cffrelpk has space allocated for deg+1 lowest doubles,
 *            for the real parts of the coefficients;
 *   cffimrtb has space allocated for deg+1 highest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimrix has space allocated for deg+1 second highest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimrmi has space allocated for deg+1 third highest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimrrg has space allocated for deg+1 fourth highest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimrpk has space allocated for deg+1 fifth highest doubles,
 *            for the imaginary parts of the coefficients.
 *   cffimltb has space allocated for deg+1 fifth lowest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimlix has space allocated for deg+1 fourth lowest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimlmi has space allocated for deg+1 third lowest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimlrg has space allocated for deg+1 second lowest doubles,
 *            for the imaginary parts of the coefficients;
 *   cffimlpk has space allocated for deg+1 lowest doubles,
 *            for the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx      nvr integers in the range from 0 to dim-1,
 *            idx(k) is the index of the k-th variable in the monomial;
 *   exp      nvr positive integers with the powers of the variables,
 *            exp(k) is the power of the variable with index idx(k);
 *   cffrertb holds deg+1 highest doubles of the real parts
 *            of the coefficients of the power series;
 *   cffrerix holds deg+1 second highest doubles of the real parts
 *            of the coefficients of the power series;
 *   cffrermi holds deg+1 third highest doubles of the real parts
 *            of the coefficients of the power series;
 *   cffrerrg holds deg+1 fourth highest doubles of the real parts
 *            of the coefficients of the power series;
 *   cffrerpk holds deg+1 fifth highest doubles of the real parts
 *            of the coefficients of the power series;
 *   cffimltb holds deg+1 fifth lowest doubles of the imaginary parts
 *            of the coefficients of the power series;
 *   cffimlix holds deg+1 fourth lowest doubles of the imaginary parts
 *            of the coefficients of the power series;
 *   cffimlmi holds deg+1 third lowest doubles of the imaginary parts
 *            of the coefficients of the power series;
 *   cffimlrg holds deg+1 second lowest doubles of the imaginary parts
 *            of the coefficients of the power series;
 *   cffimlpk holds deg+1 lowest doubles of the imaginary parts
 *            of the coefficients of the power series. */

void make_real10_input
 ( int dim, int deg,
   double **datartb, double **datarix, double **datarmi,
   double **datarrg, double **datarpk,
   double **dataltb, double **datalix, double **datalmi,
   double **datalrg, double **datalpk );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim       dimension of the input;
 *   deg       degree of the power series;
 *   datartb   space allocated for dim arrays of deg+1 doubles;
 *   datarix   space allocated for dim arrays of deg+1 doubles;
 *   datarmi   space allocated for dim arrays of deg+1 doubles;
 *   datarrg   space allocated for dim arrays of deg+1 doubles;
 *   datarpk   space allocated for dim arrays of deg+1 doubles;
 *   dataltb   space allocated for dim arrays of deg+1 doubles;
 *   datalix   space allocated for dim arrays of deg+1 doubles;
 *   datalmi   space allocated for dim arrays of deg+1 doubles;
 *   datalrg   space allocated for dim arrays of deg+1 doubles;
 *   datalpk   space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datartb  stores the highest doubles of the input series;
 *   datarix  stores the second highest doubles of the input series;
 *   datarmi  stores the third highest doubles of the input series;
 *   datarrg  stores the fourth highest doubles of the input series;
 *   datarpk  stores the fifth highest doubles of the input series;
 *   dataltb  stores the fifth lowest doubles of the input series;
 *   datalix  stores the fourth lowest doubles of the input series;
 *   datalmi  stores the third lowest doubles of the input series;
 *   datalrg  stores the second lowest doubles of the input series;
 *   datalpk  stores the lowest doubles of the input series;
 *            data[i][j] is the j-th coefficient of the i-th series
 *            for i in 0..dim-1 and j in 0..deg. */

void make_complex10_input
 ( int dim, int deg,
   double **datarertb, double **datarerix, double **datarermi,
   double **datarerrg, double **datarerpk,
   double **datareltb, double **datarelix, double **datarelmi,
   double **datarelrg, double **datarelpk,
   double **dataimrtb, double **dataimrix, double **dataimrmi,
   double **dataimrrg, double **dataimrpk,
   double **dataimltb, double **dataimlix, double **dataimlmi,
   double **dataimlrg, double **dataimlpk );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim       dimension of the input;
 *   deg       degree of the power series;
 *   datarertb has space allocated for the highest doubles of the real parts
 *             of dim series of degree deg.
 *   datarerix has space allocated for the second highest doubles of the
 *             real parts of dim series of degree deg.
 *   datarermi has space allocated for the third highest doubles of the
 *             real parts of dim series of degree deg.
 *   datarerrg has space allocated for the fourth highest doubles of the
 *             real parts of dim series of degree deg;
 *   datarerpk has space allocated for the fifth highest doubles of the
 *             real parts of dim series of degree deg;
 *   datareltb has space allocated for the fifth lowest doubles of the
 *             real parts of dim series of degree deg.
 *   datarelix has space allocated for the fourth lowest doubles of the
 *             real parts of dim series of degree deg.
 *   datarelmi has space allocated for the third lowest doubles of the
 *             real parts of dim series of degree deg.
 *   datarelrg has space allocated for the second lowest doubles of the
 *             real parts of dim series of degree deg;
 *   datarelpk has space allocated for the lowest doubles of the real parts
 *             of dim series of degree deg;
 *   dataimrtb has space allocated for the highest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimrix has space allocated for the second highest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimrmi has space allocated for the middle doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimrrg has space allocated for the second lowest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimrpk has space allocated for the lowest doubles of the
 *             imaginary parts of dim series of degree deg.
 *   dataimltb has space allocated for the fifth lowest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimlix has space allocated for the fourth lowest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimlmi has space allocated for the third lowest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimlrg has space allocated for the second lowest doubles of the
 *             imaginary parts of dim series of degree deg;
 *   dataimlpk has space allocated for the lowest doubles of the
 *             imaginary parts of dim series of degree deg.
 *
 * ON RETURN :
 *   datarertb stores the highest doubles of the real parts of the data,
 *   datarerix stores the second highest doubles of the real parts,
 *   datarermi stores the third highest doubles of the real parts,
 *   datarerrg stores the fourth highest doubles of the real parts,
 *   datarerpk stores the fifth highest doubles of the real parts,
 *   datareltb stores the fifth lowest doubles of the real parts,
 *   datarelix stores the fourth lowest doubles of the real parts,
 *   datarelmi stores the third lowest doubles of the real parts,
 *   datarelrg stores the second lowest doubles of the real parts,
 *   datarelpk stores the lowest doubles of the real parts of the data,
 *   dataimrtb stores the highest doubles of the imaginary parts of the data,
 *   dataimrix stores the second highest doubles of the imaginary parts,
 *   dataimrmi stores the third highest doubles of the imaginary parts,
 *   dataimrrg stores the fourth highest doubles of the imaginary parts,
 *   dataimrpk stores the fifth highest doubles of the imaginary parts,
 *   dataimltb stores the fifth lowest doubles of the imaginary parts,
 *   dataimlix stores the fourth lowest doubles of the imaginary parts,
 *   dataimlmi stores the third lowest doubles of the imaginary parts,
 *   dataimlrg stores the second lowest doubles of the imaginary parts,
 *   dataimlpk stores the lowest doubles of the imaginary parts of the data,
 *             data[i][j] is the j-th coefficient of the i-th series,
 *             for i in 0..dim-1 and j in 0..deg. */

#endif
