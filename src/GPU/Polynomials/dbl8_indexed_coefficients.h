// The file dbl8_indexed_coefficients.h specifies functions to define
// coefficients of indexed systems in octo double precision.

#ifndef __dbl8_indexed_coefficients_h__
#define __dbl8_indexed_coefficients_h__

int dbl8_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo, int vrblvl );
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
 *   csthihihi has space for dim pointers;
 *   cstlohihi has space for dim pointers;
 *   csthilohi has space for dim pointers;
 *   cstlolohi has space for dim pointers;
 *   csthihilo has space for dim pointers;
 *   cstlohilo has space for dim pointers;
 *   csthilolo has space for dim pointers;
 *   cstlololo has space for dim pointers;
 *   cffhihihi has space for dim pointers;
 *   cfflohihi has space for dim pointers;
 *   cffhilohi has space for dim pointers;
 *   cfflolohi has space for dim pointers;
 *   cffhihilo has space for dim pointers;
 *   cfflohilo has space for dim pointers;
 *   cffhilolo has space for dim pointers;
 *   cfflololo has space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   csthihihi are the highest doubles of the constant;
 *   cstlohihi are the second highest doubles of the constant;
 *   csthilohi are the third highest doubles of the constant;
 *   cstlohihi are the fourth highest doubles of the constant;
 *   csthihilo are the fourth lowest doubles of the constant;
 *   cstlohilo are the third lowest doubles of the constant;
 *   csthilolo are the second lowest doubles of the constant;
 *   cstlololo are the lowest doubles of the constant;
 *   cffhihihi are the highest doubles of all monomial coefficients;
 *   cfflohihi are the second highest doubles of all monomial coefficients;
 *   cffhilohi are the third highest doubles of all monomial coefficients;
 *   cffhihilo are the fourth highest doubles of all monomial coefficients;
 *   cffhihilo are the fourth lowest doubles of all monomial coefficients;
 *   cfflohilo are the third lowest doubles of all monomial coefficients;
 *   cffhilolo are the second lowest doubles of all monomial coefficients;
 *   cfflololo are the lowest doubles of all monomial coefficients. */

int cmplx8_make_coefficients
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihihi, double **cstrelohihi,
   double **cstrehilohi, double **cstrelolohi,
   double **cstrehihilo, double **cstrelohilo,
   double **cstrehilolo, double **cstrelololo,
   double **cstimhihihi, double **cstimlohihi,
   double **cstimhilohi, double **cstimlolohi,
   double **cstimhihilo, double **cstimlohilo,
   double **cstimhilolo, double **cstimlololo,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo, int vrblvl );
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
 *   cstrehihihi has space for dim pointers;
 *   cstrelohihi has space for dim pointers;
 *   cstrehilohi has space for dim pointers;
 *   cstrelolohi has space for dim pointers;
 *   cstrehihilo has space for dim pointers;
 *   cstrelohilo has space for dim pointers;
 *   cstrehilolo has space for dim pointers;
 *   cstrelololo has space for dim pointers;
 *   cstimhihihi has space for dim pointers;
 *   cstimlohihi has space for dim pointers;
 *   cstimhilohi has space for dim pointers;
 *   cstimlolohi has space for dim pointers;
 *   cstimhihilo has space for dim pointers;
 *   cstimlohilo has space for dim pointers;
 *   cstimhilolo has space for dim pointers;
 *   cstimlololo has space for dim pointers;
 *   cffrehihihi has space for dim pointers;
 *   cffrelohihi has space for dim pointers;
 *   cffrehilohi has space for dim pointers;
 *   cffrelolohi has space for dim pointers;
 *   cffrehihilo has space for dim pointers;
 *   cffrelohilo has space for dim pointers;
 *   cffrehilolo has space for dim pointers;
 *   cffrelololo has space for dim pointers;
 *   cffimhihihi has space for dim pointers;
 *   cffimlohihi has space for dim pointers;
 *   cffimhilohi has space for dim pointers;
 *   cffimlolohi has space for dim pointers;
 *   cffimhihilo has space for dim pointers;
 *   cffimlohilo has space for dim pointers;
 *   cffimhilolo has space for dim pointers;
 *   cffimlololo has space for dim pointers;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   cstrehihihi are the highest real doubles of the constant;
 *   cstrelohihi are the second highest real doubles of the constant;
 *   cstrehilohi are the third highest real doubles of the constant;
 *   cstrelolohi are the fourth highest real doubles of the constant;
 *   cstrehihilo are the fourth lowest real doubles of the constant;
 *   cstrelohilo are the third lowest real doubles of the constant;
 *   cstrehilolo are the second lowest real doubles of the constant;
 *   cstrelololo are the lowest real doubles of the constant;
 *   cstimhihihi are the highest imaginary parts of the constant;
 *   cstimlohihi are the second highest imaginary parts of of the constant;
 *   cstimhilohi are the third highest imaginary parts of of the constant;
 *   cstimlolohi are the fourth highest imaginary parts of of the constant;
 *   cstimhihilo are the fourth lowest imaginary parts of the constant;
 *   cstimlohilo are the third lowest imaginary parts of the constant;
 *   cstimhilolo are the second lowest imaginary parts of the constant;
 *   cstimlololo are the lowest imaginary parts of of the constant;
 *   cffrehihihi are the highest real parts of the coefficients;
 *   cffrelohihi are the second highest real parts of the coefficients;
 *   cffrehilohi are the third highest real parts of the coefficients;
 *   cffrelolohi are the fourth highest real parts of the coefficients;
 *   cffrehihilo are the fourth lowest real parts of the coefficients;
 *   cffrelohilo are the third lowest real parts of the coefficients;
 *   cffrehilolo are the second lowest real parts of the coefficients;
 *   cffrelololo are the lowest real parts of the coefficients;
 *   cffimhihihi are the highest imaginary parts of the coefficients;
 *   cffimlohihi are the second highest imag parts of the coefficients;
 *   cffimhilohi are the third highest imag parts of the coefficients;
 *   cffimlolohi are the fourth highest imag parts of the coefficients;
 *   cffimhihilo are the fourth lowest imag parts of the coefficients;
 *   cffimlohilo are the third lowest imag parts of the coefficients;
 *   cffimhilolo are the second lowest imag parts of the coefficients;
 *   cffimlololo are the lowest imag parts of the coefficients. */

#endif
