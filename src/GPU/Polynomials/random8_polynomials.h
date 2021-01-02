// The file random8_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in octo double precision.

#ifndef __random8_polynomials_h__
#define __random8_polynomials_h__

bool make_real8_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo );
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

void make_real8_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo );
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

void make_real8_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg.
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
#endif
