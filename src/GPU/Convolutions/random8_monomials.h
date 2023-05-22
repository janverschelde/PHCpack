// The file random8_monomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a monomial in octo double precision.

#ifndef __random8_monomials_h__
#define __random8_monomials_h__

bool make_real8_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nvr       number of variables with positive power in the monomial;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   exp       space allocated for nvr integers;
 *   cffhihihi has space allocated for deg+1 doubles;
 *   cfflohihi has space allocated for deg+1 doubles;
 *   cffhilohi has space allocated for deg+1 doubles;
 *   cfflolohi has space allocated for deg+1 doubles;
 *   cffhihilo has space allocated for deg+1 doubles;
 *   cfflohilo has space allocated for deg+1 doubles;
 *   cffhilolo has space allocated for deg+1 doubles;
 *   cfflololo has space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx       nvr integers in the range from 0 to dim-1,
 *             idx(k) is the index of the k-th variable in the monomial;
 *   exp       nvr positive integers with the powers of the variables,
 *             exp(k) is the power of the variable with index idx(k);
 *   cffhihihi stores the deg+1 highest coefficient doubles;
 *   cfflohihi stores the deg+1 second highest coefficient doubles;
 *   cffhilohi stores the deg+1 third highest coefficient doubles;
 *   cfflolohi stores the deg+1 fourth highest coefficient doubles;
 *   cffhihilo stores the deg+1 fourth lowest coefficient doubles;
 *   cfflohilo stores the deg+1 third highest coefficient doubles;
 *   cffhilolo stores the deg+1 second lowest coefficient doubles;
 *   cfflololo stores the deg+1 lowest coefficient doubles. */

void random_octo_complex
 ( double *rehihihi, double *relohihi, double *rehilohi, double *relolohi,
   double *rehihilo, double *relohilo, double *rehilolo, double *relololo,
   double *imhihihi, double *imlohihi, double *imhilohi, double *imlolohi,
   double *imhihilo, double *imlohilo, double *imhilolo, double *imlololo );
/*
 * DESCRIPTION :
 *   Generates a random complex number on the unit circle.
 *
 * ON RETURN :
 *   rehihihi   the highest double of the real part;
 *   relohihi   the second highest double of the real part;
 *   rehilohi   the third highest double of the real part;
 *   relolohi   the fourth highest double of the real part;
 *   rehihilo   the fourth lowest double of the real part;
 *   relohilo   the third lowest double of the real part;
 *   rehilolo   the second lowest double of the real part;
 *   relololo   the lowest double of the real part;
 *   imhihihi   the highest double of the imag part;
 *   imlohihi   the second highest double of the imag part;
 *   imhilohi   the third highest double of the imag part;
 *   imlolohi   the fourth highest double of the imag part;
 *   imhihilo   the fourth lowest double of the imag part;
 *   imlohilo   the third lowest double of the imag part;
 *   imhilolo   the second lowest double of the imag part;
 *   imlololo   the lowest double of the imag part. */

bool make_complex8_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim         dimension, total number of variables;
 *   nvr         number of variables with positive power in the monomial;
 *   pwr         largest power of a variable;
 *   deg         degree of the power series coefficient;
 *   exp         space allocated for nvr integers;
 *   cffrehihihi has space allocated for deg+1 doubles, for the
 *               highest doubles of the real parts of the coefficients;
 *   cffrelohihi has space allocated for deg+1 doubles, for the second
 *               highest doubles of the real parts of the coefficients;
 *   cffrehilohi has space allocated for deg+1 doubles, for the third
 *               highest doubles of the real parts of the coefficients;
 *   cffrelolohi has space allocated for deg+1 doubles, for the fourth
 *               highest doubles of the real parts of the coefficients;
 *   cffrehihilo has space allocated for deg+1 doubles, for the fourth
 *               lowest doubles of the real parts of the coefficients;
 *   cffrelohilo has space allocated for deg+1 doubles, for the third
 *               lowest doubles of the real parts of the coefficients;
 *   cffrehilolo has space allocated for deg+1 doubles, for the second
 *               lowest doubles of the real parts of the coefficients;
 *   cffrelololo has space allocated for deg+1 doubles, for the lowest
 *               doubles of the real parts of the coefficients;
 *   cffimhihihi has space allocated for deg+1 doubles, for the highest
 *               doubles of the imaginary parts of the coefficients;
 *   cffimlohihi has space allocated for deg+1 doubles, for the second
 *               highest doubles of the imaginary parts of the coefficients;
 *   cffimhilohi has space allocated for deg+1 doubles, for the third
 *               highest doubles of the imaginary parts of the coefficients;
 *   cffimlolohi has space allocated for deg+1 doubles, for the fourth
 *               highest doubles of the imaginary parts of the coefficients;
 *   cffimhihilo has space allocated for deg+1 doubles, for the fourth
 *               lowest doubles of the imaginary parts of the coefficients.
 *   cffimlohilo has space allocated for deg+1 doubles, for the third
 *               lowest doubles of the imaginary parts of the coefficients.
 *   cffimhilolo has space allocated for deg+1 doubles, for the second
 *               lowest doubles of the imaginary parts of the coefficients.
 *   cffimlololo has space allocated for deg+1 doubles, for the lowest
 *               doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx         nvr integers in the range from 0 to dim-1,
 *               idx(k) is the index of the k-th variable in the monomial;
 *   exp         nvr positive integers with the powers of the variables,
 *               exp(k) is the power of the variable with index idx(k);
 *   cffrehihihi holds the deg+1 highest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrelohihi holds the deg+1 second highest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrehilohi holds the deg+1 third highest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrelolohi holds the deg+1 fourth highest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrehihilo holds the deg+1 fourth lowest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrelohilo holds the deg+1 third lowest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrehilolo holds the deg+1 second lowest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffrelololo holds the deg+1 lowest doubles of the real parts
 *               of the coefficients of the power series;
 *   cffimhihihi holds the deg+1 highest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimlohihi holds the deg+1 second highest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimhilohi holds the deg+1 third highest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimlolohi holds the deg+1 fourth highest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimhihilo holds the deg+1 fourth lowest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimlohilo holds the deg+1 third lowest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimhilolo holds the deg+1 second lowest doubles of the imaginary parts
 *               of the coefficients of the power series;
 *   cffimlololo holds deg+1 lowest doubles of the imaginary parts
 *               of the coefficients of the power series. */

void make_real8_input
 ( int dim, int deg,
   double **datahihihi, double **datalohihi,
   double **datahilohi, double **datalolohi,
   double **datahihilo, double **datalohilo,
   double **datahilolo, double **datalololo );
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
 *   datahihihi has space allocated for dim arrays of deg+1 doubles;
 *   datalohihi has space allocated for dim arrays of deg+1 doubles;
 *   datahilohi has space allocated for dim arrays of deg+1 doubles;
 *   datalolohi has space allocated for dim arrays of deg+1 doubles.
 *   datahihilo has space allocated for dim arrays of deg+1 doubles;
 *   datalohilo has space allocated for dim arrays of deg+1 doubles;
 *   datahilolo has space allocated for dim arrays of deg+1 doubles;
 *   datalololo has space allocated for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   datahihihi holds the highest doubles of the input series,
 *   datalohihi holds the second highest doubles of the input series,
 *   datahilohi holds the third highest doubles of the input series,
 *   datalolohi holds the fourth highest doubles of the input series,
 *   datahihilo holds the fourth lowest doubles of the input series,
 *   datalohilo holds the third lowest doubles of the input series,
 *   datahilolo holds the second lowest doubles of the input series,
 *   datalololo holds the lowest doubles of the input series,
 *              data[i][j] is the j-th coefficient of the i-th series,
 *              for i in 0..dim-1 and j in 0..deg;
 *              datalohihi[i][j] is the highest double of the j-th 
 *              coefficient of the i-th series,
 *              for i in 0..dim-1 and j in 0..deg. */

void make_complex8_input
 ( int dim, int deg,
   double **datarehihihi, double **datarelohihi,
   double **datarehilohi, double **datarelolohi,
   double **datarehihilo, double **datarelohilo,
   double **datarehilolo, double **datarelololo,
   double **dataimhihihi, double **dataimlohihi,
   double **dataimhilohi, double **dataimlolohi,
   double **dataimhihilo, double **dataimlohilo,
   double **dataimhilolo, double **dataimlololo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *   The input series are such that the odd indexed series
 *   are the inverses of the previous even indexed series.
 *   In case of an odd dimension, the last series equals one.
 *   The complete product of all series must thus equal one.
 *
 * ON ENTRY :
 *   dim          dimension of the input;
 *   deg          degree of the power series;
 *   datarehihihi has space allocated for the highest doubles
 *                of the real parts of dim series of degree deg.
 *   datarelohihi has space allocated for the second highest doubles
 *                of the real parts of dim series of degree deg.
 *   datarehilohi has space allocated for the third highest doubles
 *                of the real parts of dim series of degree deg.
 *   datarelolohi has space allocated for the fourth highest doubles
 *                of the real parts of dim series of degree deg.
 *   datarehihilo has space allocated for the fourth lowest doubles
 *                of the real parts of dim series of degree deg;
 *   datarelohilo has space allocated for the third lowest doubles
 *                of the real parts of dim series of degree deg;
 *   datarehilolo has space allocated for the second lowest doubles
 *                of the real parts of dim series of degree deg;
 *   datarelololo has space allocated for the lowest doubles
 *                of the real parts of dim series of degree deg;
 *   dataimhihihi has space allocated for the highest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimlohihi has space allocated for the second highest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimhilohi has space allocated for the third highest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimlolohi has space allocated for the fourth highest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimhihilo has space allocated for the fourth lowest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimlohilo has space allocated for the third lowest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimhilolo has space allocated for the second lowest doubles
 *                of the imaginary parts of dim series of degree deg;
 *   dataimlololo has space allocated for the lowest doubles
 *                of the imaginary parts of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehihihi stores the highest doubles of the real parts,
 *   datarelohihi stores the second highest doubles of the real parts,
 *   datarehilohi stores the third highest doubles of the real parts,
 *   datarelolohi stores the fourth highest doubles of the real parts,
 *   datarehihilo stores the fourth lowest doubles of the real parts,
 *   datarelohilo stores the third lowest doubles of the real parts,
 *   datarehilolo stores the second lowest doubles of the real parts,
 *   datarelololo stores the lowest doubles of the real parts,
 *   dataimhihihi stores the highest doubles of the imaginary parts,
 *   dataimlohihi stores the second highest doubles of the imaginary parts,
 *   dataimhilohi stores the third highest doubles of the imaginary parts,
 *   dataimlolohi stores the fourth highest doubles of the imaginary parts,
 *   dataimhihilo stores the fourth lowest doubles of the imaginary parts,
 *   dataimlohilo stores the third lowest doubles of the imaginary parts,
 *   dataimhilolo stores the second lowest doubles of the imaginary parts,
 *   dataimlololo stores the lowest doubles of the imaginary parts,
 *                data[i][j] is the j-th coefficient of the i-th series,
 *                for i in 0..dim-1 and j in 0..deg. */

#endif
