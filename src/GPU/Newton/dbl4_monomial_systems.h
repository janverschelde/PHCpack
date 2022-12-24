// The file dbl4_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in quad double precision.

#ifndef __dbl4_monomial_systems_h__
#define __dbl4_monomial_systems_h__

void make_real4_exponentials
 ( int dim, int  deg,
   double **shihi, double **slohi, double **shilo, double **slolo );
/*
 * DESCRIPTION :
 *   Returns the expansions of exp(c*x) for random coefficients c,
 *   for c in the union of the intervals [-2, -1] and [1, 2].
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   shihi    space for dim arrays of size deg+1;
 *   slohi    space for dim arrays of size deg+1;
 *   shilo    space for dim arrays of size deg+1;
 *   slolo    space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   shihi    highest doubles of the series expansions;
 *   slohi    second highest doubles of the series expansions;
 *   shilo    second lowest doubles of the series expansions;
 *   slolo    lowest doubles of the series expansions. */

void make_complex4_exponentials
 ( int dim, int deg,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo );
/*
 * DESCRIPTION :
 *   Returns dim expansions of exp(x) truncated at degree deg,
 *   for random complex values of x on the complex unit circle.
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   srehihi  space for dim arrays of size deg+1;
 *   srelohi  space for dim arrays of size deg+1;
 *   srehilo  space for dim arrays of size deg+1;
 *   srelolo  space for dim arrays of size deg+1;
 *   simhihi  space for dim arrays of size deg+1;
 *   simlohi  space for dim arrays of size deg+1;
 *   simhilo  space for dim arrays of size deg+1;
 *   simlolo  space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   srehihi  highest doubles of the real parts of the series;
 *   srelohi  second highest doubles of the real parts of the series;
 *   srehilo  second lowest doubles of the real parts of the series;
 *   srelolo  lowest doubles of the real parts of the series;
 *   simhihi  highest doubles of the imaginary parts of the series;
 *   simlohi  second highest doubles of the imaginary parts of the series;
 *   simhilo  second lowest doubles of the imaginary parts of the series;
 *   simlolo  lowest doubles of the imaginary parts of the series. */

void evaluate_real4_monomials
 ( int dim, int deg, int **rowsA,
   double **shihi, double **slohi, double **shilo, double **slolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   shihi    highest doubles of the series;
 *   slohi    second highest doubles of the series;
 *   shilo    second lowest doubles of the series;
 *   slolo    lowest doubles of the series;
 *   rhshihi  has space for dim arrays of size deg+1;
 *   rhslohi  has space for dim arrays of size deg+1;
 *   rhshilo  has space for dim arrays of size deg+1;
 *   rhslolo  has space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhshihi  are the highest doubles of the real parts of the evaluations;
 *   rhslohi  are the 2nd highest doubles of the real parts of rhs;
 *   rhshilo  are the 2nd lowest doubles of the real parts of rhs;
 *   rhslolo  are the lowest doubles of the real parts of rhs. */

void evaluate_complex4_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihi, double **srelohi, double **srehilo, double **srelolo,
   double **simhihi, double **simlohi, double **simhilo, double **simlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   srehihi  highest doubles of the real parts of the series;
 *   srelohi  second highest doubles of the real parts of the series;
 *   srehilo  second lowest doubles of the real parts of the series;
 *   srelolo  lowest doubles of the real parts of the series;
 *   simhihi  highest doubles of the imaginary parts of the series;
 *   simlohi  second highest doubles of the imaginary parts of the series;
 *   simhilo  second lowest doubles of the imaginary parts of the series;
 *   simlolo  lowest doubles of the imaginary parts of the series;
 *   rhsrehihi has space for dim arrays of size deg+1;
 *   rhsrelohi has space for dim arrays of size deg+1;
 *   rhsrehilo has space for dim arrays of size deg+1;
 *   rhsrelolo has space for dim arrays of size deg+1;
 *   rhsimhihi has space for dim arrays of size deg+1;
 *   rhsimlohi has space for dim arrays of size deg+1;
 *   rhsimhilo has space for dim arrays of size deg+1;
 *   rhsimlolo has space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhsrehihi are the highest doubles of the real parts of the evaluations;
 *   rhsrelohi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelolo are the lowest doubles of the real parts of rhs;
 *   rhsimhihi are the highest doubles of the imaginary parts of the rhs;
 *   rhsimlohi are the 2nd highest doubles of the imaginary parts of the rhs;
 *   rhsimhilo are the 2nd lowest doubles of the imaginary parts of the rhs;
 *   rhsimlolo are the lowest doubles of the imaginary parts of the rhs. */

void make_real4_coefficients
 ( int nbrcol, int dim,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo );
/*
 * DESCRIPTION :
 *   Generates random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cffhihi  space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflohi  space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffhilo  space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cfflolo  space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cffhihi  highest doubles of the coefficients;
 *   cfflohi  second highest doubles of the coefficients;
 *   cffhilo  second lowest doubles of the coefficients;
 *   cfflolo  lowest doubles of the coefficients,
 *   cff      cff[i] has the coefficients for the i-th column,
 *            cff[i][j] is the coefficient of the j-th monomial
 *            in the i-th column. */

void make_complex4_coefficients
 ( int nbrcol, int dim,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo );
/*
 * DESCRIPTION :
 *   Generates complex random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cffrehihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrehilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffrelolo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhihi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlohi has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimhilo has space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffimlolo has space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cffrehihi are the highest doubles of the real parts of the coefficients;
 *   cffrelohi are 2nd highest doubles of the real parts of the coefficients;
 *   cffrehilo are 2nd lowest doubles of the real parts of the coefficients;
 *   cffrelolo are the lowest doubles of the real parts of the coefficients;
 *   cffre    doubles of the real parts of the coefficients,
 *            cff[i] has the real coefficient part for the i-th column,
 *            cff[i][j] is the real part of the coefficient
 *            of the j-th monomial in the i-th column;
 *   cffimhihi are the highest doubles of the imag parts of the coefficients;
 *   cffimlohi are 2nd highest doubles of the imag parts of the coefficients;
 *   cffimhilo are 2nd lowest doubles of the imag parts of the coefficients;
 *   cffimlolo are the lowest doubles of the imag parts of the coefficients,
 *   cffim    low doubles of the imaginary parts of the coefficients,
 *            cff[i] has the imag coefficient part for the i-th column,
 *            cff[i][j] is the imaginary part of the coefficient
 *            of the j-th monomial in the i-th column. */

void evaluate_real4_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   int vrblvl );
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
 *   cffhihi  highest double coefficients of the system;
 *   cfflohi  second highest double coefficients of the system;
 *   cffhilo  second lowest double coefficients of the system;
 *   cfflolo  lowest double coefficients of the system;
 *   xhihi    highest double coefficients of the series;
 *   xlohi    second highest double coefficients of the series;
 *   xhilo    second lowest double coefficients of the series;
 *   xlolo    lowest double coefficients of the series;
 *   rhshihi  space for dim arrays of size deg+1;
 *   rhslohi  space for dim arrays of size deg+1;
 *   rhshilo  space for dim arrays of size deg+1;
 *   rhslolo  space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhshihi  are the highest doubles of the real parts of the evaluations;
 *   rhslohi  are the 2nd highest doubles of the real parts of rhs;
 *   rhshilo  are the 2nd lowest doubles of the real parts of rhs;
 *   rhslolo  are the lowest doubles of the real parts of rhs. */

void evaluate_complex4_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **xrehihi, double **xrelohi, double **xrehilo, double **xrelolo,
   double **ximhihi, double **ximlohi, double **ximhilo, double **ximlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo, int vrblvl );
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
 *   cffrehihi are the highest double real parts of the coeffs of the system;
 *   cffrelohi are 2nd highest double real parts of the coeffs of the system;
 *   cffrehilo are 2nd lowest double real parts of the coeffs of the system;
 *   cffrelolo are the lowest double real parts of the coeffs of the system;
 *   cffimhihi are highest double imag parts of the coeffs of the system;
 *   cffimlohi are 2nd highest double imag parts of the coeffs of the system;
 *   cffimhilo are 2nd lowest double imag parts of the coeffs of the system;
 *   cffimlolo are the lowest double imag parts of the coeffs of the system;
 *   xrehihi  highest double real parts of the coefficients of the series;
 *   xrelohi  2nd highest double real parts of the coefficients of the series;
 *   xrehilo  2nd lowest double real parts of the coefficients of the series;
 *   xrelolo  lowest double real parts of the coefficients of the series;
 *   ximhihi  highest double imag parts of the coefficients of the series;
 *   ximlohi  2nd highest double imag parts of the coefficients of the series;
 *   ximhilo  2nd lowest double imag parts of the coefficients of the series;
 *   ximlolo  lowest double imag parts of the coefficients of the series;
 *   rhsrehihi has space for dim arrays of size deg+1;
 *   rhsrelohi has space for dim arrays of size deg+1;
 *   rhsrehilo has space for dim arrays of size deg+1;
 *   rhsrelolo has space for dim arrays of size deg+1;
 *   rhsimhihi has space for dim arrays of size deg+1;
 *   rhsimlohi has space for dim arrays of size deg+1;
 *   rhsimhilo has space for dim arrays of size deg+1;
 *   rhsimlolo has space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhsrehihi are the highest doubles of the real parts of the evaluations;
 *   rhsrelohi are the 2nd highest doubles of the real parts of rhs;
 *   rhsrehilo are the 2nd lowest doubles of the real parts of rhs;
 *   rhsrelolo are the lowest doubles of the real parts of rhs;
 *   rhsimhihi are the highest doubles of the imaginary parts of the rhs;
 *   rhsimlohi are the 2nd highest doubles of the imaginary parts of the rhs;
 *   rhsimhilo are the 2nd lowest doubles of the imaginary parts of the rhs;
 *   rhsimlolo are the lowest doubles of the imaginary parts of the rhs. */

#endif
