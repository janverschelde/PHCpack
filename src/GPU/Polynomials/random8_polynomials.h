// The file random8_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in octo double precision.

#ifndef __random8_polynomials_h__
#define __random8_polynomials_h__

bool make_real8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers with counts for the number of variables,
 *             nvr[k] has the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   csthihihi has space allocated for deg+1 doubles;
 *   cstlohihi has space allocated for deg+1 doubles;
 *   csthilohi has space allocated for deg+1 doubles;
 *   cstlolohi has space allocated for deg+1 doubles;
 *   csthihilo has space allocated for deg+1 doubles;
 *   cstlohilo has space allocated for deg+1 doubles;
 *   csthilolo has space allocated for deg+1 doubles;
 *   cstlololo has space allocated for deg+1 doubles;
 *   cffhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   csthihihi stores the deg+1 highest coefficient doubles
 *             of the constant monomial;
 *   cstlohihi stores the deg+1 second highest coefficient doubles
 *             of the constant monomial;
 *   csthilohi stores the deg+1 third highest coefficient doubles
 *             of the constant monomial;
 *   cstlolohi stores the deg+1 fourth highest coefficient doubles
 *             of the constant monomial;
 *   csthihilo stores the deg+1 fourth lowest coefficient doubles
 *             of the constant monomial;
 *   cstlohilo stores the deg+1 third highest coefficient doubles
 *             of the constant monomial;
 *   csthilolo stores the deg+1 second lowest coefficient doubles
 *             of the constant monomial;
 *   cstlololo stores the deg+1 lowest coefficient doubles
 *             of the constant monomial;
 *   cffhihihi stores the highest parts of the series for each monomial;
 *   cfflohihi stores the second highest parts of the series;
 *   cffhilolo stores the third highest parts of the series;
 *   cfflolohi stores the fourth highest parts of the series;
 *   cfflohihi stores the fourth lowest parts of the series;
 *   cfflolohi stores the third lowest parts of the series;
 *   cfflohilo stores the second lowest parts of the series;
 *   cfflololo stores the lowest parts of the series for each monomial. */

bool make_cmplx8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers with counts for the number of variables,
 *             nvr[k] has the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   cstrehihihi has space allocated for deg+1 doubles;
 *   cstrelohihi has space allocated for deg+1 doubles;
 *   cstrehilohi has space allocated for deg+1 doubles;
 *   cstrelolohi has space allocated for deg+1 doubles;
 *   cstrehihilo has space allocated for deg+1 doubles;
 *   cstrelohilo has space allocated for deg+1 doubles;
 *   cstrehilolo has space allocated for deg+1 doubles;
 *   cstrelololo has space allocated for deg+1 doubles;
 *   cstimhihihi has space allocated for deg+1 doubles;
 *   cstimlohihi has space allocated for deg+1 doubles;
 *   cstimhilohi has space allocated for deg+1 doubles;
 *   cstimlolohi has space allocated for deg+1 doubles;
 *   cstimhihilo has space allocated for deg+1 doubles;
 *   cstimlohilo has space allocated for deg+1 doubles;
 *   cstimhilolo has space allocated for deg+1 doubles;
 *   cstimlololo has space allocated for deg+1 doubles;
 *   cffrehihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelololo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   cstrehihihi stores the deg+1 highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohihi stores the deg+1 second highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilohi stores the deg+1 third highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelolohi stores the deg+1 fourth highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehihilo stores the deg+1 fourth lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohilo stores the deg+1 third lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilolo stores the deg+1 second lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelololo stores the deg+1 lowest real coefficient doubles
 *             of the constant monomial;
 *   cstimhihihi stores the deg+1 highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohihi stores the deg+1 second highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilohi stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlolohi stores the deg+1 fourth highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhihilo stores the deg+1 fourth lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohilo stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilolo stores the deg+1 second lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlololo stores the deg+1 lowest imag coefficient doubles
 *             of the constant monomial;
 *   cffrehihihi stores the highest real parts of the series for each monomial;
 *   cffrehilohi stores the second highest real parts of the series;
 *   cffrehihilo stores the third highest real parts of the series;
 *   cffrehilolo stores the fourth highest real parts of the series;
 *   cffrelohihi stores the fourth lowest real parts of the series;
 *   cffrelolohi stores the third lowest real parts of the series;
 *   cffrelohilo stores the second lowest real parts of the series;
 *   cffrelololo stores the lowest real parts of the series;
 *   cffimhihihi stores the highest imag parts of the series;
 *   cffimhilohi stores the second highest imag parts of the series;
 *   cffimhihilo stores the third highest imag parts of the series;
 *   cffimhilolo stores the fourth highest imag parts of the series;
 *   cffimlohihi stores the fourth lowest imag parts of the series;
 *   cffimlolohi stores the third lowest imag parts of the series;
 *   cffimlohilo stores the second lowest imag parts of the series;
 *   cffimlololo stores the lowest imag parts of the series. */

void make_real8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nbr index vectors;
 *   csthihihi has space allocated for deg+1 doubles;
 *   cstlohihi has space allocated for deg+1 doubles;
 *   csthilohi has space allocated for deg+1 doubles;
 *   cstlolohi has space allocated for deg+1 doubles;
 *   csthihilo has space allocated for deg+1 doubles;
 *   cstlohilo has space allocated for deg+1 doubles;
 *   csthilolo has space allocated for deg+1 doubles;
 *   cstlololo has space allocated for deg+1 doubles;
 *   cffhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   csthihihi stores the deg+1 highest coefficient doubles
 *             of the constant monomial;
 *   cstlohihi stores the deg+1 second highest coefficient doubles
 *             of the constant monomial;
 *   csthilohi stores the deg+1 third highest coefficient doubles
 *             of the constant monomial;
 *   cstlolohi stores the deg+1 fourth highest coefficient doubles
 *             of the constant monomial;
 *   csthihilo stores the deg+1 fourth lowest coefficient doubles
 *             of the constant monomial;
 *   cstlohilo stores the deg+1 third highest coefficient doubles
 *             of the constant monomial;
 *   csthilolo stores the deg+1 second lowest coefficient doubles
 *             of the constant monomial;
 *   cstlololo stores the deg+1 lowest coefficient doubles
 *             of the constant monomial;
 *   cffhihihi stores the highest parts of the series for each monomial,
 *             cffhihihi[k] has the deg+1 highest coefficients of monomial k;
 *   cffhilohi stores the second highest parts of the series, cffhilohi[k]
 *             has the deg+1 second highest coefficients of monomial k;
 *   cffhihilo stores the third highest parts of the series, cffhihilo[k]
 *             has the deg+1 third highest coefficients of monomial k;
 *   cffhilolo stores the fourth highest parts of the series, cffhilolo[k]
 *             has the deg+1 fourth highest coefficients of monomial k;
 *   cfflohihi stores the fourth lowest parts of the series, cfflohihi[k]
 *             has the deg+1 fourth lowest coefficients of monomial k;
 *   cfflolohi stores the third lowest parts of the series, cfflolohi[k]
 *             has the deg+1 third lowest coefficients of monomial k;
 *   cfflohilo stores the second lowest parts of the series, cfflohilo[k]
 *             has the deg+1 second lowest coefficients of monomial k;
 *   cfflololo stores the lowest parts of the series for each monomial,
 *             cfflololo[k] has the deg+1 lowest coefficients of monomial k. */

void make_cmplx8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg,
 *   for complex data.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nva index vectors;
 *   cstrehihihi has space allocated for deg+1 doubles;
 *   cstrelohihi has space allocated for deg+1 doubles;
 *   cstrehilohi has space allocated for deg+1 doubles;
 *   cstrelolohi has space allocated for deg+1 doubles;
 *   cstrehihilo has space allocated for deg+1 doubles;
 *   cstrelohilo has space allocated for deg+1 doubles;
 *   cstrehilolo has space allocated for deg+1 doubles;
 *   cstrelololo has space allocated for deg+1 doubles;
 *   cstimhihihi has space allocated for deg+1 doubles;
 *   cstimlohihi has space allocated for deg+1 doubles;
 *   cstimhilohi has space allocated for deg+1 doubles;
 *   cstimlolohi has space allocated for deg+1 doubles;
 *   cstimhihilo has space allocated for deg+1 doubles;
 *   cstimlohilo has space allocated for deg+1 doubles;
 *   cstimhilolo has space allocated for deg+1 doubles;
 *   cstimlololo has space allocated for deg+1 doubles;
 *   cffrehihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelololo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   cstrehihihi stores the deg+1 highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohihi stores the deg+1 second highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilohi stores the deg+1 third highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelolohi stores the deg+1 fourth highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehihilo stores the deg+1 fourth lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohilo stores the deg+1 third lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilolo stores the deg+1 second lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelololo stores the deg+1 lowest real coefficient doubles
 *             of the constant monomial;
 *   cstimhihihi stores the deg+1 highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohihi stores the deg+1 second highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilohi stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlolohi stores the deg+1 fourth highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhihilo stores the deg+1 fourth lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohilo stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilolo stores the deg+1 second lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlololo stores the deg+1 lowest imag coefficient doubles
 *             of the constant monomial;
 *   cffrehihihi stores the highest real parts of the series for each monomial;
 *   cffrehilohi stores the second highest real parts of the series;
 *   cffrehihilo stores the third highest real parts of the series;
 *   cffrehilolo stores the fourth highest real parts of the series;
 *   cffrelohihi stores the fourth lowest real parts of the series;
 *   cffrelolohi stores the third lowest real parts of the series;
 *   cffrelohilo stores the second lowest real parts of the series;
 *   cffrelololo stores the lowest real parts of the series;
 *   cffimhihihi stores the highest imag parts of the series;
 *   cffimhilohi stores the second highest imag parts of the series;
 *   cffimhihilo stores the third highest imag parts of the series;
 *   cffimhilolo stores the fourth highest imag parts of the series;
 *   cffimlohihi stores the fourth lowest imag parts of the series;
 *   cffimlolohi stores the third lowest imag parts of the series;
 *   cffimlohilo stores the second lowest imag parts of the series;
 *   cffimlololo stores the lowest imag parts of the series. */

void make_real8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg,
 *   for real data.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nva index vectors;
 *   csthihihi has space allocated for deg+1 doubles;
 *   cstlohihi has space allocated for deg+1 doubles;
 *   csthilohi has space allocated for deg+1 doubles;
 *   cstlolohi has space allocated for deg+1 doubles;
 *   csthihilo has space allocated for deg+1 doubles;
 *   cstlohilo has space allocated for deg+1 doubles;
 *   csthilolo has space allocated for deg+1 doubles;
 *   cstlololo has space allocated for deg+1 doubles;
 *   cffhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cfflohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cfflololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   csthihihi stores the deg+1 highest coefficient doubles
 *             of the constant monomial;
 *   cstlohihi stores the deg+1 second highest coefficient doubles
 *             of the constant monomial;
 *   csthilohi stores the deg+1 third highest coefficient doubles
 *             of the constant monomial;
 *   cstlolohi stores the deg+1 fourth highest coefficient doubles
 *             of the constant monomial;
 *   csthihilo stores the deg+1 fourth lowest coefficient doubles
 *             of the constant monomial;
 *   cstlohilo stores the deg+1 third lowest coefficient doubles
 *             of the constant monomial;
 *   csthilolo stores the deg+1 second lowest coefficient doubles
 *             of the constant monomial;
 *   cstlololo stores the deg+1 lowest coefficient doubles
 *             of the constant monomial;
 *   cffhihihi stores the highest parts of the series for each monomial,
 *             cffhihihi[k] has the deg+1 highest coefficients of monomial k;
 *   cffhilohi stores the second highest parts of the series, cffhilohi[k]
 *             has the deg+1 second highest coefficients of monomial k;
 *   cffhihilo stores the third highest parts of the series, cffhihilo[k]
 *             has the deg+1 third highest coefficients of monomial k;
 *   cffhilolo stores the fourth highest parts of the series, cffhilolo[k]
 *             has the deg+1 fourth highest coefficients of monomial k;
 *   cfflohihi stores the fourth lowest parts of the series, cfflohihi[k]
 *             has the deg+1 fourth lowest coefficients of monomial k;
 *   cfflolohi stores the third lowest parts of the series, cfflolohi[k]
 *             has the deg+1 third lowest coefficients of monomial k;
 *   cfflohilo stores the second lowest parts of the series, cfflohilo[k]
 *             has the deg+1 second lowest coefficients of monomial k;
 *   cfflololo stores the lowest parts of the series for each monomial,
 *             cfflololo[k] has the deg+1 lowest coefficients of monomial k. */

void make_cmplx8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg,
 *   for complex data.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nva index vectors;
 *   cstrehihihi has space allocated for deg+1 doubles;
 *   cstrelohihi has space allocated for deg+1 doubles;
 *   cstrehilohi has space allocated for deg+1 doubles;
 *   cstrelolohi has space allocated for deg+1 doubles;
 *   cstrehihilo has space allocated for deg+1 doubles;
 *   cstrelohilo has space allocated for deg+1 doubles;
 *   cstrehilolo has space allocated for deg+1 doubles;
 *   cstrelololo has space allocated for deg+1 doubles;
 *   cstimhihihi has space allocated for deg+1 doubles;
 *   cstimlohihi has space allocated for deg+1 doubles;
 *   cstimhilohi has space allocated for deg+1 doubles;
 *   cstimlolohi has space allocated for deg+1 doubles;
 *   cstimhihilo has space allocated for deg+1 doubles;
 *   cstimlohilo has space allocated for deg+1 doubles;
 *   cstimhilolo has space allocated for deg+1 doubles;
 *   cstimlololo has space allocated for deg+1 doubles;
 *   cffrehihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelololo has space allocated for nbr arrays of deg+1 doubles.
 *   cffimhihihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlolohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlololo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   cstrehihihi stores the deg+1 highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohihi stores the deg+1 second highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilohi stores the deg+1 third highest real coefficient doubles
 *             of the constant monomial;
 *   cstrelolohi stores the deg+1 fourth highest real coefficient doubles
 *             of the constant monomial;
 *   cstrehihilo stores the deg+1 fourth lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelohilo stores the deg+1 third lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrehilolo stores the deg+1 second lowest real coefficient doubles
 *             of the constant monomial;
 *   cstrelololo stores the deg+1 lowest real coefficient doubles
 *             of the constant monomial;
 *   cstimhihihi stores the deg+1 highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohihi stores the deg+1 second highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilohi stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlolohi stores the deg+1 fourth highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhihilo stores the deg+1 fourth lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlohilo stores the deg+1 third highest imag coefficient doubles
 *             of the constant monomial;
 *   cstimhilolo stores the deg+1 second lowest imag coefficient doubles
 *             of the constant monomial;
 *   cstimlololo stores the deg+1 lowest imag coefficient doubles
 *             of the constant monomial;
 *   cffrehihihi stores the highest real parts of the series for each monomial;
 *   cffrehilohi stores the second highest real parts of the series;
 *   cffrehihilo stores the third highest real parts of the series;
 *   cffrehilolo stores the fourth highest real parts of the series;
 *   cffrelohihi stores the fourth lowest real parts of the series;
 *   cffrelolohi stores the third lowest real parts of the series;
 *   cffrelohilo stores the second lowest real parts of the series;
 *   cffrelololo stores the lowest real parts of the series;
 *   cffimhihihi stores the highest imag parts of the series;
 *   cffimhilohi stores the second highest imag parts of the series;
 *   cffimhihilo stores the third highest imag parts of the series;
 *   cffimhilolo stores the fourth highest imag parts of the series;
 *   cffimlohihi stores the fourth lowest imag parts of the series;
 *   cffimlolohi stores the third lowest imag parts of the series;
 *   cffimlohilo stores the second lowest imag parts of the series;
 *   cffimlololo stores the lowest imag parts of the series. */

#endif
