// The file dbl8_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in octo double precision.

#ifndef __dbl8_monomial_systems_h__
#define __dbl8_monomial_systems_h__

void make_real8_exponentials
 ( int dim, int  deg,
   double **shihihi, double **slohihi, double **shilohi, double **slolohi,
   double **shihilo, double **slohilo, double **shilolo, double **slololo );
/*
 * DESCRIPTION :
 *   Returns the expansions of exp(c*x) for random coefficients c,
 *   for c in the union of the intervals [-2, -1] and [1, 2].
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   shihihi  space for dim arrays of size deg+1;
 *   slohihi  space for dim arrays of size deg+1;
 *   shilohi  space for dim arrays of size deg+1;
 *   slolohi  space for dim arrays of size deg+1;
 *   shihilo  space for dim arrays of size deg+1;
 *   slohilo  space for dim arrays of size deg+1;
 *   shilolo  space for dim arrays of size deg+1;
 *   slololo  space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   shihihi  highest doubles of the series expansions;
 *   slohihi  second highest doubles of the series expansions;
 *   shilohi  third highest doubles of the series expansions;
 *   slolohi  fourth highest doubles of the series expansions;
 *   shihilo  fourth lowest doubles of the series expansions;
 *   slohilo  third lowest doubles of the series expansions;
 *   shilolo  second lowest doubles of the series expansions;
 *   slololo  lowest doubles of the series expansions. */

void make_complex8_exponentials
 ( int dim, int deg,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo );
/*
 * DESCRIPTION :
 *   Returns dim expansions of exp(x) truncated at degree deg,
 *   for random complex values of x on the complex unit circle.
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   srehihihi has space for dim arrays of size deg+1;
 *   srelohihi has space for dim arrays of size deg+1;
 *   srehihihi has space for dim arrays of size deg+1;
 *   srelohihi has space for dim arrays of size deg+1;
 *   srehilolo has space for dim arrays of size deg+1;
 *   srelololo has space for dim arrays of size deg+1;
 *   srehilolo has space for dim arrays of size deg+1;
 *   srelololo has space for dim arrays of size deg+1;
 *   simhihihi has space for dim arrays of size deg+1;
 *   simlohihi has space for dim arrays of size deg+1;
 *   simhilohi has space for dim arrays of size deg+1;
 *   simlolohi has space for dim arrays of size deg+1.
 *   simhihilo has space for dim arrays of size deg+1;
 *   simlohilo has space for dim arrays of size deg+1;
 *   simhilolo has space for dim arrays of size deg+1;
 *   simlololo has space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   srehihihi are the highest doubles of the real parts of the series;
 *   srelohihi are the 2nd highest doubles of the real parts of the series;
 *   srehilohi are the 3rd highest doubles of the real parts of the series;
 *   srelolohi are the 4rth highest doubles of the real parts of the series;
 *   srehihilo are the 4rth lowest doubles of the real parts of the series;
 *   srelohilo are the 3rd lowest doubles of the real parts of the series;
 *   srehilolo are the 2nd lowest doubles of the real parts of the series;
 *   srelololo are the lowest doubles of the real parts of the series;
 *   simhihihi are the highest doubles of the imag parts of the series;
 *   simlohihi are the 2nd highest doubles of the imag parts of the series;
 *   simhilohi are the 3rd highest doubles of the imag parts of the series;
 *   simlolohi are the 4th highest doubles of the imag parts of the series;
 *   simhihilo are the 4th lowest doubles of the imag parts of the series;
 *   simlohilo are the 3rd lowest doubles of the imag parts of the series;
 *   simhilolo are the 2nd lowest doubles of the imag parts of the series;
 *   simlololo are the lowest doubles of the imag parts of the series. */

void evaluate_real8_monomials
 ( int dim, int deg, int **rowsA,
   double **shihihi, double **slohihi, double **shilohi, double **slolohi,
   double **shihilo, double **slohilo, double **shilolo, double **slololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   shihihi  are the highest doubles of the real parts of the series;
 *   slohihi  are the 2nd highest doubles of the series;
 *   shilohi  are the 3rd highest doubles of the series;
 *   slolohi  are the 4th highest doubles of the series;
 *   shihilo  are the 4th lowest doubles of the series;
 *   slohilo  are the 3rd lowest doubles of the series;
 *   shilolo  are the 2nd lowest doubles of the series;
 *   slololo  are the lowest doubles of the series;
 *   rhshihihi has space for dim arrays of size deg+1;
 *   rhslohihi has space for dim arrays of size deg+1;
 *   rhshilohi has space for dim arrays of size deg+1;
 *   rhslolohi has space for dim arrays of size deg+1;
 *   rhshihilo has space for dim arrays of size deg+1;
 *   rhslohilo has space for dim arrays of size deg+1;
 *   rhshilolo has space for dim arrays of size deg+1;
 *   rhslololo has space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhshihihi are the highest doubles of the real parts of the evaluations;
 *   rhslohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhshilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhslolohi are the 4th highest doubles of the real parts of rhs;
 *   rhshihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhslohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhshilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhslololo are the lowest doubles of the real parts of rhs. */

void evaluate_complex8_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   srehihihi are the highest doubles of the real parts of the series;
 *   srelohihi are the 2nd highest doubles of the real parts of the series;
 *   srehilohi are the 3rd highest doubles of the real parts of the series;
 *   srelolohi are the 4th highest doubles of the real parts of the series;
 *   srehihilo are the 4th lowest doubles of the real parts of the series;
 *   srelohilo are the 3rd lowest doubles of the real parts of the series;
 *   srehilolo are the 2nd lowest doubles of the real parts of the series;
 *   srelololo are the lowest doubles of the real parts of the series;
 *   simhihihi are the highest doubles of the imag parts of the series;
 *   simlohihi are the 2nd highest doubles of the imag parts of the series;
 *   simhilohi are the 3rd highest doubles of the imag parts of the series;
 *   simlolohi are the 4th highest doubles of the imag parts of the series;
 *   simhihilo are the 4th lowest doubles of the imag parts of the series;
 *   simlohilo are the 3rd lowest doubles of the imag parts of the series;
 *   simhilolo are the 2nd lowest doubles of the imag parts of the series;
 *   simlololo are the lowest doubles of the imag parts of the series;
 *   rhsrehihihi has space for dim arrays of size deg+1;
 *   rhsrelohihi has space for dim arrays of size deg+1;
 *   rhsrehilohi has space for dim arrays of size deg+1;
 *   rhsrelolohi has space for dim arrays of size deg+1;
 *   rhsrehihilo has space for dim arrays of size deg+1;
 *   rhsrelohilo has space for dim arrays of size deg+1;
 *   rhsrehilolo has space for dim arrays of size deg+1;
 *   rhsrelololo has space for dim arrays of size deg+1;
 *   rhsimhihihi has space for dim arrays of size deg+1;
 *   rhsimlohihi has space for dim arrays of size deg+1;
 *   rhsimhilohi has space for dim arrays of size deg+1;
 *   rhsimlolohi has space for dim arrays of size deg+1.
 *   rhsimhihilo has space for dim arrays of size deg+1;
 *   rhsimlohilo has space for dim arrays of size deg+1;
 *   rhsimhilolo has space for dim arrays of size deg+1;
 *   rhsimlololo has space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhsrehihihi are the highest doubles of the real parts of the evaluations;
 *   rhsrelohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhsrelolohi are the 4th highest doubles of the real parts of rhs;
 *   rhsrehihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhsrelohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhsrehilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelololo are the lowest doubles of the real parts of rhs;
 *   rhsimhihihi are the highest doubles of the imag parts of the rhs;
 *   rhsimlohihi are the 2nd highest doubles of the imag parts of the rhs;
 *   rhsimhilohi are the 3rd highest doubles of the imag parts of the rhs;
 *   rhsimlolohi are the 4th highest doubles of the imag parts of the rhs;
 *   rhsimhihilo are the 4th lowest doubles of the imag parts of the rhs;
 *   rhsimlohilo are the 3rd lowest doubles of the imag parts of the rhs;
 *   rhsimhilolo are the 2nd lowest doubles of the imag parts of the rhs;
 *   rhsimlololo are the lowest doubles of the imag parts of the rhs. */

void make_real8_coefficients
 ( int nbrcol, int dim,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo );
/*
 * DESCRIPTION :
 *   Generates random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cffhihihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflohihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffhilohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflolohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffhihilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflohilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffhilolo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflololo has space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cffhihihi are the highest doubles of the coefficients;
 *   cfflohihi are the second highest doubles of the coefficients;
 *   cfflohihi are the third highest doubles of the coefficients;
 *   cfflolohi are the fourth highest doubles of the coefficients;
 *   cffhihilo are the fourth lowest doubles of the coefficients;
 *   cfflohilo are the third lowest doubles of the coefficients,
 *   cffhilolo are the second lowest doubles of the coefficients;
 *   cfflololo are the lowest doubles of the coefficients,
 *   cff      cff[i] has the coefficients for the i-th column,
 *            cff[i][j] is the coefficient of the j-th monomial
 *            in the i-th column. */

void make_complex8_coefficients
 ( int nbrcol, int dim,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo );
/*
 * DESCRIPTION :
 *   Generates complex random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cffrehihihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelohihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrehilohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelolohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrehihilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelohilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrehilolo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelololo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhihihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlohihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhilohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlolohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhihilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlohilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhilolo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlololo has space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cffrehihihi are the highest doubles of the real parts of the coeffs;
 *   cffrelohihi are 2nd highest doubles of the real parts of the coeffs;
 *   cffrehilohi are 3rd highest doubles of the real parts of the coeffs;
 *   cffrelolohi are 4th highest doubles of the real parts of the coeffs;
 *   cffrehihilo are 4th lowest doubles of the real parts of the coeffs;
 *   cffrelohilo are 3rd lowest doubles of the real parts of the coeffs;
 *   cffrehilolo are 2nd lowest doubles of the real parts of the coeffs;
 *   cffrelololo are the lowest doubles of the real parts of the coeffs;
 *   cffre    doubles of the real parts of the coefficients,
 *            cff[i] has the real coefficient part for the i-th column,
 *            cff[i][j] is the real part of the coefficient
 *            of the j-th monomial in the i-th column;
 *   cffimhihihi are the highest doubles of the imag parts of the coeffs;
 *   cffimlohihi are 2nd highest doubles of the imag parts of the coeffs;
 *   cffimhilohi are 3rd highest doubles of the imag parts of the coeffs;
 *   cffimlolohi are 4th highest doubles of the imag parts of the coeffs,
 *   cffimhihilo are 4th lowest doubles of the imag parts of the coeffs;
 *   cffimlohilo are 3rd lowest doubles of the imag parts of the coeffs;
 *   cffimhilolo are 2nd lowest doubles of the imag parts of the coeffs;
 *   cffimlololo are the lowest doubles of the imag parts of the coeffs,
 *   cffim    low doubles of the imaginary parts of the coefficients,
 *            cff[i] has the imag coefficient part for the i-th column,
 *            cff[i][j] is the imaginary part of the coefficient
 *            of the j-th monomial in the i-th column. */

void evaluate_real8_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the polynomials defined by the column representation
 *   at real series to make the right hand side of a system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   nbrcol   is the number of columns;
 *   nvr      nvr[[i][j] is the number of variables of the j-th monomial
 *            in the i-th column.
 *   idx      idx[i][j][k] is the index of the k-th variable which appears
 *            in the j-th monomial of the i-th column;
 *   rowsA    matrix of dimension dim where the rows of A 
 *            are used as work space during the evaluation;
 *   cffhihihi are the highest double coefficients of the system;
 *   cfflohihi are the second highest double coefficients of the system;
 *   cffhilohi are the third highest double coefficients of the system;
 *   cfflolohi are the fourth highest double coefficients of the system;
 *   cffhihilo are the fourth lowest double coefficients of the system;
 *   cfflohilo are the third lowest double coefficients of the system;
 *   cffhilolo are the second lowest double coefficients of the system;
 *   cfflololo are the lowest double coefficients of the system;
 *   xhihihi  highest double coefficients of the series;
 *   xlohihi  second highest double coefficients of the series;
 *   xhilohi  third highest double coefficients of the series;
 *   xlolohi  fourth highest double coefficients of the series;
 *   xhihilo  fourth lowest double coefficients of the series;
 *   xlohilo  third lowest double coefficients of the series;
 *   xhilolo  second lowest double coefficients of the series;
 *   xlololo  lowest double coefficients of the series;
 *   rhshihihi has space for dim arrays of size deg+1;
 *   rhslohihi has space for dim arrays of size deg+1;
 *   rhshilohi has space for dim arrays of size deg+1;
 *   rhslolohi has space for dim arrays of size deg+1;
 *   rhshihilo has space for dim arrays of size deg+1;
 *   rhslohilo has space for dim arrays of size deg+1;
 *   rhshilolo has space for dim arrays of size deg+1;
 *   rhslololo has space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhshihihi are the highest doubles of the real parts of the evaluations;
 *   rhslohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhshilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhslolohi are the 4th highest doubles of the real parts of rhs;
 *   rhshihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhslohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhshilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhslololo are the lowest doubles of the real parts of rhs. */

void evaluate_complex8_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **xrehihihi, double **xrelohihi,
   double **xrehilohi, double **xrelolohi,
   double **xrehihilo, double **xrelohilo,
   double **xrehilolo, double **xrelololo,
   double **ximhihihi, double **ximlohihi,
   double **ximhilohi, double **ximlolohi,
   double **ximhihilo, double **ximlohilo,
   double **ximhilolo, double **ximlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the polynomials defined by the column representation
 *   at real series to make the right hand side of a system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   nbrcol   is the number of columns;
 *   nvr      nvr[[i][j] is the number of variables of the j-th monomial
 *            in the i-th column.
 *   idx      idx[i][j][k] is the index of the k-th variable which appears
 *            in the j-th monomial of the i-th column;
 *   rowsA    matrix of dimension dim where the rows of A 
 *            are used as work space during the evaluation;
 *   cffrehihihi has the highest double real parts of the sys coeffs;
 *   cffrelohihi has the 2nd highest double real parts of the sys coeffs;
 *   cffrehilohi has the 3rd highest double real parts of the sys coeffs;
 *   cffrelolohi has the 4th highest double real parts of the sys coeffs;
 *   cffrehihilo has the 4th lowest double real parts of the sys coeffs;
 *   cffrelohilo has the 3rd lowest double real parts of the sys coeffs;
 *   cffrehilolo has the 2nd lowest double real parts of the sys coeffs;
 *   cffrelololo has the lowest double real parts of the sys coeffs;
 *   cffimhihihi has the highest double imag parts of the sys coeffs;
 *   cffimlohihi has the 2nd highest double imag parts of the sys coeffs;
 *   cffimhilohi has the 3rd highest double imag parts of the sys coeffs;
 *   cffimlolohi has the 4th highest double imag parts of the sys coeffs;
 *   cffimhihilo has the 4th lowest double imag parts of the sys coeffs;
 *   cffimlohilo has the 3rd lowest double imag parts of the sys coeffs;
 *   cffimhilolo has the 2nd lowest double imag parts of the sys coeffs;
 *   cffimlololo has the lowest double imag parts of the sys coeffs;
 *   xrehihihi has the highest double real parts of the coefficients;
 *   xrelohihi has the 2nd highest double real parts of the coefficients;
 *   xrehilohi has the 3rd highest double real parts of the coefficients;
 *   xrelolohi has the 4th highest double real parts of the coefficients;
 *   xrehihilo has the 4th lowest double real parts of the coefficients;
 *   xrelohilo has the 3rd lowest double real parts of the coefficients;
 *   xrehilolo has the 2nd lowest double real parts of the coefficients;
 *   xrelololo has the lowest double real parts of the coefficients;
 *   ximhihihi has the highest double imag parts of the coefficients;
 *   ximlohihi has the 2nd highest double imag parts of the coefficients;
 *   ximhilohi has the 3rd highest double imag parts of the coefficients;
 *   ximlolohi has the 4th highest double imag parts of the coefficients;
 *   ximhihilo has the 4th lowest double imag parts of the coefficients;
 *   ximlohilo has the 3rd lowest double imag parts of the coefficients;
 *   ximhilolo has the 2nd lowest double imag parts of the coefficients;
 *   ximlololo has the lowest double imag parts of the coefficients;
 *   rhsrehihihi has space for dim arrays of size deg+1;
 *   rhsrelohihi has space for dim arrays of size deg+1;
 *   rhsrehilohi has space for dim arrays of size deg+1;
 *   rhsrelolohi has space for dim arrays of size deg+1;
 *   rhsrehihilo has space for dim arrays of size deg+1;
 *   rhsrelohilo has space for dim arrays of size deg+1;
 *   rhsrehilolo has space for dim arrays of size deg+1;
 *   rhsrelololo has space for dim arrays of size deg+1;
 *   rhsimhihihi has space for dim arrays of size deg+1;
 *   rhsimlohihi has space for dim arrays of size deg+1;
 *   rhsimhilohi has space for dim arrays of size deg+1;
 *   rhsimlolohi has space for dim arrays of size deg+1;
 *   rhsimhihilo has space for dim arrays of size deg+1;
 *   rhsimlohilo has space for dim arrays of size deg+1;
 *   rhsimhilolo has space for dim arrays of size deg+1;
 *   rhsimlololo has space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhsrehihihi are the highest doubles of the real parts of the evaluations;
 *   rhsrelohihi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilohi are the 3rd highest doubles of the real parts of rhs;
 *   rhsrelolohi are the 4th highest doubles of the real parts of rhs;
 *   rhsrehihilo are the 4th lowest doubles of the real parts of rhs;
 *   rhsrelohilo are the 3rd lowest doubles of the real parts of rhs;
 *   rhsrehilolo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelololo are the lowest doubles of the real parts of rhs;
 *   rhsimhihihi are the highest doubles of the imag parts of the rhs;
 *   rhsimlohihi are the 2nd highest doubles of the imag parts of the rhs;
 *   rhsimhilohi are the 3rd highest doubles of the imag parts of the rhs;
 *   rhsimlolohi are the 4th highest doubles of the imag parts of the rhs;
 *   rhsimhihilo are the 4th lowest doubles of the imag parts of the rhs;
 *   rhsimlohilo are the 3rd lowest doubles of the imag parts of the rhs;
 *   rhsimhilolo are the 2nd lowest doubles of the imag parts of the rhs;
 *   rhsimlololo are the lowest doubles of the imag parts of the rhs. */

#endif
