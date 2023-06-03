// The file dbl_indexed_coefficients.h specifies functions to define
// coefficients of indexed systems in double precision.

#ifndef __dbl_indexed_coefficients_h__
#define __dbl_indexed_coefficients_h__

int dbl_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, int vrblvl );
/*
 * DESCRIPTION :
 *   Makes the real constant and the real coefficients 
 *   of an indexed polynomial system.
 *
 * ON ENTRY :
 *   dim      number of equations and variables in the system;
 *   deg      degree of truncation of the series in the system;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial;
 *   cst      space for dim pointers;
 *   cff      space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cst      coefficients of the constant series;
 *   cff      coefficients of all monomials in the system. */

int cmplx_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstre, double **cstim, double ***cffre, double ***cffim,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Makes the complex constant and the complex coefficients
 *   of an indexed polynomial system.
 *
 * ON ENTRY :
 *   dim      number of equations and variables in the system;
 *   deg      degree of truncation of the series in the system;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial;
 *   cstre    space for dim pointers;
 *   cstim    space for dim pointers;
 *   cffre    space for dim pointers;
 *   cffim    space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cstre    real parts of the coefficients of the constant series;
 *   cstim    imaginary parts of the coefficients of the constant series;
 *   cffre    real parts of the coefficients of all monomials in the system;
 *   cffim    imaginary parts of the coefficients of all monomials. */

#endif
