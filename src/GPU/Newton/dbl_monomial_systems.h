// The file dbl_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in double precision.

#ifndef __dbl_monomial_systems_h__
#define __dbl_monomial_systems_h__

void make_real_exponentials ( int dim, int  deg, double **s );
/*
 * DESCRIPTION :
 *   Returns the expansions of exp(c*x) for random coefficients c,
 *   for c in the union of the intervals [-2, -1] and [1, 2].
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   s        space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   s         coefficients of the power series expansions. */

void make_complex_exponentials
 ( int dim, int deg, double *angles, double **sre, double **sim );
/*
 * DESCRIPTION :
 *   Returns dim expansions of exp(c*x) truncated at degree deg,
 *   for random complex values of c on the complex unit circle,
 *   where the i-th c is the defined as cos(angle[i]) + I*sin(angle[i]).
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   angles   has space for dim angles;
 *   sre      space for dim arrays of size deg+1;
 *   sim      space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   angles   contains the angles of the random complex numbers,
 *            cos(angles[i]) is sre[i][1] and sin(angles[i]) is sim[i][1];
 *   sre      real parts of the coefficients of the dim power series;
 *   sim      imaginary parts of the coefficients of the dim power series. */

void evaluate_real_monomials
 ( int dim, int deg, int **rowsA, double **x, double **rhs );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   real series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   x        coefficients of the series;
 *   rhs      space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhs      the evaluated monomials. */

void evaluate_complex_monomials
 ( int dim, int deg, int **rowsA,
   double **xre, double **xim, double **rhsre, double **rhsim );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   xre      real parts of the series;
 *   xim      imaginary parts of the series;
 *   rhsre    space for dim arrays of size deg+1;
 *   rhsim    space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhsre    real parts of the evaluated monomials;
 *   rhsim    imaginary parts of the evaluated monomials. */

void make_real_coefficients ( int nbrcol, int dim, double ***cff );
/*
 * DESCRIPTION :
 *   Generates random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cff      space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cff      cff[i] has the coefficients for the i-th column,
 *            cff[i][j] is the coefficient of the j-th monomial
 *            in the i-th column. */

void make_complex_coefficients
 ( int nbrcol, int dim, double ***cffre, double ***cffim );
/*
 * DESCRIPTION :
 *   Generates complex random double coefficients for a column system.
 *   Assigns only the leading coefficient of each coefficient series.
 *
 * ON ENTRY :
 *   nbrcol   number of columns is the leading dimension;
 *   dim      number of equations, the dimension of the system;
 *   cffre    space for nbrcol columns with at least dim doubles
 *            in each column;
 *   cffim    space for nbrcol columns with at least dim doubles
 *            in each column.
 *
 * ON RETURN :
 *   cffre    cffre[i] has the real coefficient part for the i-th column,
 *            cffre[i][j] is the real part of the coefficient
 *            of the j-th monomial in the i-th column;
 *   cffim    cffim[i] has the imag coefficient part for the i-th column,
 *            cffim[i][j] is the imaginary part of the coefficient
 *            of the j-th monomial in the i-th column. */

void evaluate_real_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cff, double **x, double **rhs, int vrblvl );
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
 *   cff      coefficients of the column system,
 *            cff[i] has the coefficients for the i-th column;
 *   x        coefficients of the series;
 *   rhs      space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhs      the evaluated monomials. */

void evaluate_complex_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffre, double ***cffim, double **xre, double **xim,
   double **rhsre, double **rhsim, int vrblvl );
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
 *   cffre    real parts of the coefficients of the column system,
 *            cffre[i] has the coefficients for the i-th column;
 *   cffim    imaginary parts of the coefficients of the column system,
 *            cffim[i] has the coefficients for the i-th column;
 *   xre      real parts of the coefficients of the series;
 *   xim      imaginary parts of the coefficients of the series;
 *   rhsre    space for dim arrays of size deg+1;
 *   rhsim    space for dim arrays of size deg+1;
 *   vrblvl   is the verbose level, if > 1, then exponents are written.
 *
 * ON RETURN :
 *   rhsre    real parts of the evaluated monomials;
 *   rhsim    imaginary parts of the evaluated monomials. */

#endif
