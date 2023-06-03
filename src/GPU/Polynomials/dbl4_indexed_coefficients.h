// The file dbl4_indexed_coefficients.h specifies functions to define
// coefficients of indexed systems in quad double precision.

#ifndef __dbl4_indexed_coefficients_h__
#define __dbl4_indexed_coefficients_h__

int dbl4_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi, double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
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
 *   csthihi  space for dim pointers;
 *   cstlohi  space for dim pointers;
 *   csthilo  space for dim pointers;
 *   cstlolo  space for dim pointers;
 *   cffhihi  space for dim pointers;
 *   cfflohi  space for dim pointers;
 *   cffhilo  space for dim pointers;
 *   cfflolo  space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   csthihi  highest doubles of the constant series;
 *   cstlohi  second highest doubles of the constant series;
 *   csthilo  second lowest doubles of the constant series;
 *   cstlolo  lowest doubles of the constant series;
 *   cffhihi  highest doubles of all monomials in the system;
 *   cfflohi  second highest doubles of all monomials in the system;
 *   cffhilo  second lowest doubles of all monomials in the system;
 *   cfflolo  lowest doubles of all monomials in the system. */

int cmplx4_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi,
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo, int vrblvl );
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
 *   cstrehihi has space for dim pointers;
 *   cstrelohi has space for dim pointers;
 *   cstrehilo has space for dim pointers;
 *   cstrelolo has space for dim pointers;
 *   cstimhihi has space for dim pointers;
 *   cstimlohi has space for dim pointers;
 *   cstimhilo has space for dim pointers;
 *   cstimlolo has space for dim pointers;
 *   cffrehihi has space for dim pointers;
 *   cffrelohi has space for dim pointers;
 *   cffrehilo has space for dim pointers;
 *   cffrelolo has space for dim pointers;
 *   cffimhihi has space for dim pointers;
 *   cffimlohi has space for dim pointers;
 *   cffimhilo has space for dim pointers;
 *   cffimlolo has space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cstrehihi are the highest real doubles of the constant;
 *   cstrelohi are the second highest real doubles of the constant;
 *   cstrehilo are the second lowest real doubles of the constant;
 *   cstrelolo are the lowest real doubles of the constant;
 *   cstimhihi are the highest imaginary parts of the constant;
 *   cstimlohi are the second highest imaginary parts of of the constant;
 *   cstimhilo are the second lowest imaginary parts of the constant;
 *   cstimlolo are the lowest imaginary parts of of the constant;
 *   cffrehihi are the highest real parts of the coefficients;
 *   cffrelohi are the second highest real parts of the coefficients;
 *   cffrehilo are the second lowest real parts of the coefficients;
 *   cffrelolo are the lowest real parts of the coefficients;
 *   cffimhihi are the highest imaginary parts of the coefficients;
 *   cffimlohi are the second highest imaginary parts of the coefficients;
 *   cffimhilo are the second lowest imaginary parts of the coefficients;
 *   cffimlolo are the lowest imaginary parts of the coefficients. */

#endif
