// The file dbl2_indexed_coefficients.h specifies functions to define
// coefficients of indexed systems in double double precision.

#ifndef __dbl2_indexed_coefficients_h__
#define __dbl2_indexed_coefficients_h__

int dbl2_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo, double ***cffhi, double ***cfflo,
   int vrblvl );
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
 *   csthi    space for dim pointers;
 *   cstlo    space for dim pointers;
 *   cffhi    space for dim pointers;
 *   cfflo    space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   csthi    high doubles of the constant series;
 *   cstlo    low doubles of the constant series;
 *   cffhi    high doubles of all monomials in the system;
 *   cfflo    low doubles of all monomials in the system. */

int cmplx2_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
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
 *   cstrehi  space for dim pointers;
 *   cstrelo  space for dim pointers;
 *   cstimhi  space for dim pointers;
 *   cstimlo  space for dim pointers;
 *   cffrehi  space for dim pointers;
 *   cffrelo  space for dim pointers;
 *   cffimhi  space for dim pointers;
 *   cffimlo  space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cstrehi  high real doubles of the coefficients of the constant series;
 *   cstrelo  low real doubles of the coefficients of the constant series;
 *   cstimhi  high imaginary parts of the coefficients of the constant series;
 *   cstimlo  low imaginary parts of the coefficients of the constant series;
 *   cffrehi  high real parts of the coefficients of all monomials;
 *   cffrelo  low real parts of the coefficients of all monomials;
 *   cffimhi  high imaginary parts of the coefficients of all monomial;
 *   cffimlo  low imaginary parts of the coefficients of all monomials. */

#endif
