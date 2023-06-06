// The file dbl8_systems_host.h specifies functions to evaluate and
// differentiate monomials with common factors in octo double precision.

#ifndef __dbl8_systems_host_h__
#define __dbl8_systems_host_h__

void CPU_dbl8_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffhihihi are the highest doubles of the monomial coefficients;
 *   cfflohihi are the 2nd highest doubles of the monomial coefficients;
 *   cffhilohi are the 3rd highest doubles of the monomial coefficients;
 *   cfflolohi are the 4th highest doubles of the monomial coefficients;
 *   cffhihilo are the 4th lowest doubles of the monomial coefficients;
 *   cfflohilo are the 3rd lowest doubles of the monomial coefficients;
 *   cffhilolo are the 2nd lowest doubles of the monomial coefficients;
 *   cfflololo are the lowest doubles of the monomial coefficients;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclohihi has space to accumulate one power series of degree deg;
 *   acchilohi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclohilo has space to accumulate one power series of degree deg;
 *   acchilolo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   inputhihihi are the highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohihi are the second highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilohi are the third highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlolohi are the fourth highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhihilo are the fourth lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohilo are the third lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilolo are the second lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlololo are the lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   outputhihihi has space for the highest output doubles;
 *   outputlohihi has space for the second highest output doubles;
 *   outputhilohi has space for the third highest output doubles;
 *   outputlolohi has space for the fourth highest output doubles;
 *   outputhihilo has space for the fourth lowest output doubles;
 *   outputlohilo has space for the third lowest output doubles;
 *   outputhilolo has space for the second lowest output doubles;
 *   outputlololo has space for the lowest output doubles;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffhihihi are the highest doubles of the evaluated factors;
 *   cfflohihi are the second highest doubles of the evaluated factors;
 *   cffhilohi are the third highest doubles of the evaluated factors;
 *   cfflolohi are the fourth highest doubles of the evaluated factors;
 *   cffhihilo are the fourth lowest doubles of the evaluated factors;
 *   cfflohilo are the third lowest doubles of the evaluated factors;
 *   cffhilolo are the second lowest doubles of the evaluated factors;
 *   cfflololo are the lowest doubles of the evaluated factors;
 *   outputhihihi are the highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohihi are the second highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilohi are the third highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlolohi are the fourth highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhihilo are the fourth lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohilo are the third lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilolo are the second lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlololo are the lowest doubles of the evaluated
 *             and differentiated monomials,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void CPU_cmplx8_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double *accrehihihi, double *accrelohihi,
   double *accrehilohi, double *accrelolohi,
   double *accrehihilo, double *accrelohilo,
   double *accrehilolo, double *accrelololo,
   double *accimhihihi, double *accimlohihi,
   double *accimhilohi, double *accimlolohi,
   double *accimhihilo, double *accimlohilo,
   double *accimhilolo, double *accimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi, 
   double **inputimhilohi, double **inputimlolohi, 
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo, 
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi, 
   double ***outputrehihilo, double ***outputrelohilo, 
   double ***outputrehilolo, double ***outputrelololo, 
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffrehihihi are the highest doubles of real parts of coefficients;
 *   cffrelohihi are the 2nd highest doubles of real parts of coefficients;
 *   cffrehilohi are the 3rd highest doubles of real parts of coefficients;
 *   cffrelolohi are the 4th highest doubles of real parts of coefficients;
 *   cffrehihilo are the 4th lowest doubles of real parts of coefficients;
 *   cffrelohilo are the 3rd lowest doubles of real parts of coefficients;
 *   cffrehilolo are the 2nd lowest doubles of real parts of coefficients;
 *   cffrelololo are the lowest doubles of real parts of coefficients;
 *   cffimhihihi are the highest doubles of imaginary parts of coefficients;
 *   cffimlohihi are the 2nd highest doubles of imag parts of coefficients;
 *   cffimhilohi are the 3rd highest doubles of imag parts of coefficients;
 *   cffimlolohi are the 4th highest doubles of imag parts of coefficients;
 *   cffimhihilo are the 4th lowest doubles of imag parts of coefficients;
 *   cffimlohilo are the 3rd lowest doubles of imag parts of coefficients;
 *   cffimhilolo are the 2nd lowest doubles of imag parts of coefficients;
 *   cffimlololo are the lowest doubles of imaginary parts of coefficients;
 *   accrehihihi has space to accumulate one power series of degree deg;
 *   accrelohihi has space to accumulate one power series of degree deg;
 *   accrehilohi has space to accumulate one power series of degree deg;
 *   accrelolohi has space to accumulate one power series of degree deg;
 *   accrehihilo has space to accumulate one power series of degree deg;
 *   accrelohilo has space to accumulate one power series of degree deg;
 *   accrehilolo has space to accumulate one power series of degree deg;
 *   accrelololo has space to accumulate one power series of degree deg;
 *   accimhihihi has space to accumulate one power series of degree deg;
 *   accimlohihi has space to accumulate one power series of degree deg;
 *   accimhilohi has space to accumulate one power series of degree deg;
 *   accimlolohi has space to accumulate one power series of degree deg;
 *   accimhihilo has space to accumulate one power series of degree deg;
 *   accimlohilo has space to accumulate one power series of degree deg;
 *   accimhilolo has space to accumulate one power series of degree deg;
 *   accimlololo has space to accumulate one power series of degree deg;
 *   inputrehihi are the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi are the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi are the third highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi are the fourth highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo are the fourth lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo are the third lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo are the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelolo are the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihi are the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi are the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi are the third highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi are the fourth highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo are the fourth lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo are the third lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo are the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolo are the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihihi has space for the output;
 *   outputrelohihi has space for the output;
 *   outputrehilohi has space for the output;
 *   outputrelolohi has space for the output;
 *   outputrehihilo has space for the output;
 *   outputrelohilo has space for the output;
 *   outputrehilolo has space for the output;
 *   outputrelololo has space for the output;
 *   outputimhihihi has space for the output;
 *   outputimlohihi has space for the output;
 *   outputimhilohi has space for the output;
 *   outputimlolohi has space for the output;
 *   outputimhihilo has space for the output;
 *   outputimlohilo has space for the output;
 *   outputimhilolo has space for the output;
 *   outputimlololo has space for the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffrehihihi are highest doubles of real parts of the common factors;
 *   cffrelohihi are 2nd highest doubles of real parts of the common factors;
 *   cffrehilohi are 3rd highest doubles of real parts of the common factors;
 *   cffrelolohi are 4th highest doubles of real parts of the common factors;
 *   cffrehihilo are 4th lowest doubles of real parts of the common factors;
 *   cffrelohilo are 3rd lowest doubles of real parts of the common factors;
 *   cffrehilolo are 2nd lowest doubles of real parts of the common factors;
 *   cffrelololo are lowest doubles of real parts of the common factors;
 *   cffimhihihi are highest doubles of imag parts of the common factors;
 *   cffimlohihi are 2nd highest doubles of imag parts of the common factors;
 *   cffimhilohi are 3rd highest doubles of imag parts of the common factors;
 *   cffimlolohi are 4th highest doubles of imag parts of the common factors;
 *   cffimhihilo are 4th lowest doubles of imag parts of the common factors;
 *   cffimlohilo are 3rd lowest doubles of imag parts of the common factors;
 *   cffimhilolo are 2nd lowest doubles of imag parts of the common factors;
 *   cffimlololo are lowest doubles of imag parts of the common factors;
 *   outputrehihihi are the highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehihihi[i][dim] is the
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehihihi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelohihi are the 2nd highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelohi[i][dim] is the 2nd 
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrehilohi are the 3rd highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehilohi[i][dim] is the 3rd 
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehilohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelolohi are the 4th highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehilohi[i][dim] is the 4th 
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelolohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrehihilo are the 4th lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehihilo[i][dim] is the 4th
 *             lowest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelohilo are the 3rd lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelohilo[i][dim] is the 3rd
 *             lowest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelohilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrehilolo are the 2nd lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehilolo[i][dim] is the 2nd
 *             lowest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehilolo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelololo are the lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelololo[i][dim] is the
 *             lowest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelololo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhihihi are the highest doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimhihihi[i][dim] is the
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhihihi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlohihi are the 2nd highest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimlohihi[i][dim] is the 2nd
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlohihi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhilohi are the 3rd highest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhilohi[i][dim] is the 3rd
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhilohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlolohi are the 4th highest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhi[i][dim] is the 4th 
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlolohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhihilo are the 4th lowest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhilo[i][dim] is the 4th
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhihilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlohilo are the 3rd lowest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhilo[i][dim] is the 3rd
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlohilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhilolo are the 2nd lowest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhilolo[i][dim] is the 2nd
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhilolo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlololo are the lowest doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimlololo[i][dim] is the
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlololo[i][idx[k]] is the derivative w.r.t. idx[k]. */

void CPU_dbl8_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **acchihihi, double **acclohihi,
   double **acchilohi, double **acclolohi,
   double **acchihilo, double **acclohilo,
   double **acchilolo, double **acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cffhihihi are the highest double coefficients of the monomials,
 *             the coefficients of the j-th monomial in the i-th column;
 *   cfflohihi are the 2nd highest double coefficients of the monomials;
 *   cffhilohi are the 3rd highest double coefficients of the monomials;
 *   cfflolohi are the 4th highest double coefficients of the monomials;
 *   cffhihilo are the 4th lowest double coefficients of the monomials;
 *   cfflohilo are the 3rd lowest double coefficients of the monomials;
 *   cffhilolo are the 2nd lowest double coefficients of the monomials;
 *   cfflololo are the lowest double coefficients of the monomials,;
 *   acchihihi has space to accumulate dim+1 power series of degree deg;
 *   acclohihi has space to accumulate dim+1 power series of degree deg;
 *   acchilohi has space to accumulate dim+1 power series of degree deg;
 *   acclolohi has space to accumulate dim+1 power series of degree deg;
 *   acchihilo has space to accumulate dim+1 power series of degree deg;
 *   acclohilo has space to accumulate dim+1 power series of degree deg;
 *   acchilolo has space to accumulate dim+1 power series of degree deg;
 *   acclololo has space to accumulate dim+1 power series of degree deg;
 *   inputhihihi are the highest double coefficients of the input;
 *   inputlohihi are the second highest double coefficients of the input;
 *   inputhilohi are the third highest double coefficients of the input;
 *   inputlolohi are the fourth highest double coefficients of the input;
 *   inputhihilo are the fourth lowest double coefficients of the input;
 *   inputlohilo are the third lowest double coefficients of the input;
 *   inputhilolo are the second lowest double coefficients of the input;
 *   inputlololo are the lowest double coefficients of the input;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihihi are the highest doubles of evaluations;
 *   funvallohihi are the second highest doubles of evaluations;
 *   funvalhilohi are the third highest doubles of evaluations;
 *   funvallolohi are the fourth highest doubles of evaluations;
 *   funvalhihilo are the fourth lowest doubles of evaluations;
 *   funvallohilo are the third lowest doubles of evaluations;
 *   funvalhilolo are the second lowest doubles of evaluations;
 *   funvallololo are the lowest doubles of evaluations;
 *   jacvalhihihi are the highest doubles of all derivatives;
 *   jacvallohihi are the second highest doubles of all derivatives;
 *   jacvalhilohi are the third highest doubles of all derivatives;
 *   jacvallolohi are the fourth highest doubles of all derivatives;
 *   jacvalhihilo are the fourth lowest doubles of all derivatives;
 *   jacvallohilo are the third lowest doubles of all derivatives;
 *   jacvalhilolo are the second lowest doubles of all derivatives;
 *   jacvallololo are the lowest doubles of all derivatives. */

void CPU_cmplx8_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **accrehihihi, double **accrelohihi,
   double **accrehilohi, double **accrelolohi,
   double **accrehihilo, double **accrelohilo,
   double **accrehilolo, double **accrelololo,
   double **accimhihihi, double **accimlohihi,
   double **accimhilohi, double **accimlolohi,
   double **accimhihilo, double **accimlohilo,
   double **accimhilolo, double **accimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates the monomials in the column representation of a system,
 *   at power series, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns;
 *   nvr       nvr[[i][j] is the number of variables of the j-th monomial
 *             in the i-th column;
 *   idx       idx[i][j][k] is the index of the k-th variable which appears
 *             in the j-th monomial of the i-th column;
 *   cffrehihihi are the highest double real parts of the coefficients,
 *             of the j-th mononmial in the i-th column;
 *   cffrelohihi are the 2nd highest double real parts of the coefficients;
 *   cffrehilohi are the 3rd highest double real parts of the coefficients;
 *   cffrelolohi are the 4th highest double real parts of the coefficients;
 *   cffrehihilo are the 4th lowest double real parts of the coefficients;
 *   cffrelohilo are the 3rd lowest double real parts of the coefficients;
 *   cffrehilolo are the 2nd lowest double real parts of the coefficients;
 *   cffrelololo are the lowest double real parts of the coefficients;
 *   cffimhihihi are the highest double imag parts of the coefficients;
 *   cffimlohihi are the 2nd highest double imag parts of the coefficients;
 *   cffimhilohi are the 3rd highest double imag parts of the coefficients;
 *   cffimlolohi are the 4th highest double imag parts of the coefficients;
 *   cffimhihilo are the 4th lowest double imag parts of the coefficients;
 *   cffimlohilo are the 3rd lowest double imag parts of the coefficients;
 *   cffimhilolo are the 2nd lowest double imag parts of the coefficients;
 *   cffimlololo are the lowest double imaginary parts of the coefficients;
 *   accrehihihi has space to accumulate dim+1 power series of degree deg;
 *   accrelohihi has space to accumulate dim+1 power series of degree deg;
 *   accrehilohi has space to accumulate dim+1 power series of degree deg;
 *   accrelolohi has space to accumulate dim+1 power series of degree deg;
 *   accrehihilo has space to accumulate dim+1 power series of degree deg;
 *   accrelohilo has space to accumulate dim+1 power series of degree deg;
 *   accrehilohi has space to accumulate dim+1 power series of degree deg;
 *   accrelololo has space to accumulate dim+1 power series of degree deg;
 *   accimhihihi has space to accumulate dim+1 power series of degree deg;
 *   accimlohihi has space to accumulate dim+1 power series of degree deg;
 *   accimhilohi has space to accumulate dim+1 power series of degree deg;
 *   accimlolohi has space to accumulate dim+1 power series of degree deg;
 *   accimhihilo has space to accumulate dim+1 power series of degree deg;
 *   accimlohilo has space to accumulate dim+1 power series of degree deg;
 *   accimhilohi has space to accumulate dim+1 power series of degree deg;
 *   accimlololo has space to accumulate dim+1 power series of degree deg;
 *   inputrehihihi are the highest double real parts of the input;
 *   inputrelohihi are the second highest double real parts of the input;
 *   inputrehilohi are the third highest double real parts of the input;
 *   inputrelolohi are the fourth highest double real parts of the input;
 *   inputrehihilo are the fourth lowest double real parts of the input;
 *   inputrelohilo are the third lowest double real parts of the input;
 *   inputrehilolo are the second lowest double real parts of the input;
 *   inputrelololo are the lowest double real parts of the input;
 *   inputimhihihi are the highest double imaginary parts of the input;
 *   inputimlohihi are the second highest double imag parts of the input;
 *   inputimhilohi are the third highest double imag parts of the input;
 *   inputimlolohi are the fourth highest double imag parts of the input;
 *   inputimhihilo are the fourth lowest double imag parts of the input;
 *   inputimlohilo are the third lowest double imag parts of the input;
 *   inputimhilolo are the second lowest double imag parts of the input;
 *   inputimlololo are the lowest double imaginary parts of the input;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihihi are highest doubles of real parts of the evaluations;
 *   funvalrelohihi are 2nd highest doubles of real parts of the evaluations;
 *   funvalrehilohi are 3rd highest doubles of real parts of the evaluations;
 *   funvalrelolohi are 4th highest doubles of real parts of the evaluations;
 *   funvalrehihilo are 4th lowest doubles of real parts of the evaluations;
 *   funvalrelohilo are 3rd lowest doubles of real parts of the evaluations;
 *   funvalrehilolo are 2nd lowest doubles of real parts of the evaluations;
 *   funvalrelololo are lowest doubles of real parts of the evaluations;
 *   funvalimhihihi are highest doubles of imag parts of the evaluations;
 *   funvalimlohihi are 2nd highest doubles of imag parts of the evaluations;
 *   funvalimhilohi are 3rd highest doubles of imag parts of the evaluations;
 *   funvalimlolohi are 4th highest doubles of imag parts of the evaluations;
 *   funvalimhihilo are 4th lowest doubles of imag parts of the evaluations;
 *   funvalimlohilo are 3rd lowest doubles of imag parts of the evaluations;
 *   funvalimhilolo are 2nd lowest doubles of imag parts of the evaluations;
 *   funvalimlololo are lowest doubles of imag parts of the evaluations;
 *   jacvalrehihihi are highest doubles of real parts of all derivatives;
 *   jacvalrelohihi are 2nd highest doubles of real parts of all derivatives;
 *   jacvalrehilohi are 3rd highest doubles of real parts of all derivatives;
 *   jacvalrelohihi are 4th highest doubles of real parts of all derivatives;
 *   jacvalrehihilo are 4th lowest doubles of real parts of all derivatives;
 *   jacvalrelohilo are 3rd lowest doubles of real parts of all derivatives;
 *   jacvalrehilolo are 2nd lowest doubles of real parts of all derivatives;
 *   jacvalrelololo are lowest doubles of real parts of all derivatives;
 *   jacvalimhihihi are highest doubles of imag parts of all derivatives;
 *   jacvalimlohihi are 2nd highest doubles of imag parts of all derivatives;
 *   jacvalimlohihi are 3rd highest doubles of imag parts of all derivatives;
 *   jacvalimhilohi are 4th highest doubles of imag parts of all derivatives;
 *   jacvalimhihilo are 4th lowest doubles of imag parts of all derivatives;
 *   jacvalimlohilo are 3rd lowest doubles of imag parts of all derivatives;
 *   jacvalimhilolo are 2nd lowest doubles of imag parts of all derivatives;
 *   jacvalimlololo are lowest doubles of imag parts of all derivatives. */

void dbl8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double damper,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi, 
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo, 
   double **funvalhilolo, double **funvallololo, 
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Linearizes the output of the evaluation and differentiation
 *   of the monomials.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   nvr       number of variables that occur in each monomial;
 *   idx       for each monomials the indices of each variable;
 *   mbhihihi  highest doubles of the right hand side of monomial system;
 *   mblohihi  second highest doubles of the right hand side;
 *   mbhilohi  third highest doubles of the right hand side;
 *   mblolohi  fourth highest doubles of the right hand side;
 *   mbhihilo  fourth lowest doubles of the right hand side;
 *   mblohilo  third lowest doubles of the right hand side;
 *   mbhilolo  second lowest doubles of the right hand side;
 *   mblololo  lowest doubles of the right hand side;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputhihihi are the highest doubles of the output
 *             of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputlohihi are the second highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilohi are the third highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlolohi are the fourth highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhihilo are the fourth lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlohilo are the third lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilolo are the second lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlololo are the lowest doubles
 *             of the output of the evaluation and differentiation;
 *   funvalhihihi has space allocated for dim power series;
 *   funvallohihi has space allocated for dim power series;
 *   funvalhilohi has space allocated for dim power series;
 *   funvallolohi has space allocated for dim power series;
 *   funvalhihilo has space allocated for dim power series;
 *   funvallohilo has space allocated for dim power series;
 *   funvalhilolo has space allocated for dim power series;
 *   funvallololo has space allocated for dim power series;
 *   rhshihihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslolohih has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhshihilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslololo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalhihihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallolohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhihilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallololo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihihi has the highest doubles of output[i][dim];
 *   funvallohihi has the second highest doubles of output[i][dim];
 *   funvallohihi has the third highest doubles of output[i][dim];
 *   funvallolohi has the fourth highest doubles of output[i][dim];
 *   funvalhihilo has the fourth lowest doubles of output[i][dim];
 *   funvallohilo has the third lowest doubles of output[i][dim];
 *   funvallohilo has the second lowest doubles of output[i][dim];
 *   funvallololo has the lowest doubles of output[i][dim];
 *   rhshihihi are the linearized right hand side are the highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohihi are the linearized right hand side are the 2nd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilohi are the linearized right hand side are the 3rd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslolohi are the linearized right hand side are the 4th highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshihilo are the linearized right hand side are the 4th lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohilo are the linearized right hand side are the 3th lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilolo are the linearized right hand side are the 2nd lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslololo are the linearized right hand side are the lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   jacvalhihihi is a series with matrices as coefficients, of highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvallohihi is a series with matrices as coefficients, of second highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilohi is a series with matrices as coefficients, of third highest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallolohi is a series with matrices as coefficients, of fourth highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhihilo is a series with matrices as coefficients, of fourth lowest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvallohilo is a series with matrices as coefficients, of third lowest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilolo is a series with matrices as coefficients, of second lowest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallololo is a series with matrices as coefficients, of lowest
 *             doubles, the leading coefficient is the Jacobian matrix. */

void cmplx8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo, double damper,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Linearizes the output of the evaluation and differentiation
 *   of the monomials, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   nvr       number of variables that occur in each monomial;
 *   idx       for each monomials the indices of each variable;
 *   mbrehihihi are the highest real parts of the right hand side series;
 *   mbrelohihi are the 2nd highest real parts of the right hand side series;
 *   mbrehilohi are the 3rd highest real parts of the right hand side series;
 *   mbrelolohi are the 4th highest real parts of the right hand side series;
 *   mbrehihilo are the 4th lowest real parts of the right hand side series;
 *   mbrelohilo are the 3rd lowest real parts of the right hand side series;
 *   mbrehilolo are the 2nd lowest real parts of the right hand side series;
 *   mbrelololo are the lowest real parts of the right hand side series;
 *   mbimhihihi are the highest imaginary parts of the right hand side series;
 *   mbimlohihi are the 2nd highest imag parts of the right hand side series;
 *   mbimhilohi are the 3rd highest imag parts of the right hand side series;
 *   mbimlolohi are the 4th highest imag parts of the right hand side series;
 *   mbimhihilo are the 4th lowest imag parts of the right hand side series;
 *   mbimlohilo are the 3rd lowest imag parts of the right hand side series;
 *   mbimhilolo are the 2nd lowest imag parts of the right hand side series;
 *   mbimlololo are the lowest imaginary parts of the right hand side series;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputrehihihi are the highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehihihi[i][dim] is the power series value, 
 *             outputrehihihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelohihi are the 2nd highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelohihi[i][dim] is the power series value, 
 *             outputrelohihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrehilohi are the 3rd highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehilohi[i][dim] is the power series value, 
 *             outputrehilohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelolohi are the 4th highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelolohi[i][dim] is the power series value, 
 *             outputrelolohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrehihilo are the 4th lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehihilo[i][dim] is the power series value, 
 *             outputrehihilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelohilo are the 3nd lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelohilo[i][dim] is the power series value, 
 *             outputrelohilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrehilolo are the 2nd lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehilolo[i][dim] is the power series value, 
 *             outputrehilolo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelololo are the lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelololo[i][dim] is the power series value, 
 *             outputrelololo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhihihi are the highest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhihihi[i][dim] is the power series value, 
 *             outputimhihihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlohihi are the 2nd highest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlohihi[i][dim] is the power series value, 
 *             outputimlohihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhilohi are the 3rd highest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhilohi[i][dim] is the power series value, 
 *             outputimhilohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlolohi are the 4th highest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlolohi[i][dim] is the power series value, 
 *             outputimlolohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhihilo are the 4th lowest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhihilo[i][dim] is the power series value, 
 *             outputimhihilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlohilo are the 3rd lowest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlohilo[i][dim] is the power series value, 
 *             outputimlohilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhilolo are the 2nd lowest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhilolo[i][dim] is the power series value, 
 *             outputimhilolo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlololo are the lowest doubles of the imaginary parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlololo[i][dim] is the power series value, 
 *             outputimlololo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   funvalrehihihi has space allocated for dim power series;
 *   funvalrelohihi has space allocated for dim power series;
 *   funvalrehilohi has space allocated for dim power series;
 *   funvalrelolohi has space allocated for dim power series;
 *   funvalrehihilo has space allocated for dim power series;
 *   funvalrelohilo has space allocated for dim power series;
 *   funvalrehilolo has space allocated for dim power series;
 *   funvalrelololo has space allocated for dim power series;
 *   funvalimhihihi has space allocated for dim power series;
 *   funvalimlohihi has space allocated for dim power series;
 *   funvalimhilohi has space allocated for dim power series;
 *   funvalimlolohi has space allocated for dim power series;
 *   funvalimhihilo has space allocated for dim power series;
 *   funvalimlohilo has space allocated for dim power series;
 *   funvalimhilolo has space allocated for dim power series;
 *   funvalimlololo has space allocated for dim power series;
 *   rhsrehihihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrehilohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelolohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrehihilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrehilolo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelololo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhihihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhilohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlolohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhihilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhilolo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlololo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalrehihihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelohihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrehilohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelolohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrehihilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelohilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrehilolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelololo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhihihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlohihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhilohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlolohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhihilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlohilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhilolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlololo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihihi are the highest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrelohihi are the second highest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrehilohi are the third highest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrelolohi are the fourth highest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrehihilo are the fourth lowest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrelohilo are the third lowest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrehilolo are the second lowest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrelololo are the lowest doubles of the real parts
 *             of the output[i][dim];
 *   funvalimhihihi are the highest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimlohihi are the second highest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimhilohi are the third highest doubles of the imaginary parts
 *             of the output[i][dim];
 *   funvalimlolohi are the fourth highest doubles of the imaginary parts
 *             of the output[i][dim];
 *   funvalimhihilo are the fourth lowest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimlohilo are the third lowest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimhilolo are the second lowest doubles of the imaginary parts
 *             of the output[i][dim];
 *   funvalimlololo are the lowest doubles of the imaginary parts
 *             of the output[i][dim];
 *   rhsrehihihi are the highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelohihi are the second highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrehilohi are the third highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelolohi are the fourth highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrehihilo are the fourth lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelohilo are the third lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrehilolo are the second lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelololo are the lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsimhihihi are the highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlohihi are the second highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimhilohi are the third highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlolohi are the fourth highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimhihilo are the fourth lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlohilo are the third lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimhilolo are the second lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlololo are the lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   jacvalrehihihi are the highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelohihi are the second highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrehilohi are the third highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelolohi are the fourth highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrehihilo are the fourth lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelohilo are the third lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalrehilolo are the second lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelololo are the lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalimhihihi are the highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimlohihi are the second highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimhilohi are the third highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimlolohi are the fourth highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimhihilo are the fourth lowest doubles imaginary parts
 *             of the matrix series;
 *   jacvalimlohilo are the third lowest doubles imaginary parts
 *             of the matrix series;
 *   jacvalimhilolo are the second lowest doubles imaginary parts
 *             of the matrix series;
 *   jacvalimlololo are the lowest doubles imaginary parts
 *             of the matrix series. */

void dbl8_define_rhs
 ( int dim, int degp1,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the right hand side for the test with the column representation
 *   of the cyclic n-roots system, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   mbhihihi  are the highest doubles of the right hand side series;
 *   mblohihi  are the second highest doubles of the right hand side series;
 *   mbhilohi  are the third highest doubles of the right hand side series;
 *   mblolohi  are the fourth highest doubles of the right hand side series;
 *   mbhihilo  are the fourth lowest doubles of the right hand side series;
 *   mblohilo  are the third lowest doubles of the right hand side series;
 *   mbhilolo  are the second lowest doubles of the right hand side series;
 *   mblololo  are the lowest doubles of the right hand side series;
 *   funvalhihihi are the highest doubles of the evaluated series;
 *   funvallohihi are the second highest doubles of the evaluated series;
 *   funvalhilohi are the third highest doubles of the evaluated series;
 *   funvallolohi are the fourth highest doubles of the evaluated series;
 *   funvalhihilo are the fourth lowest doubles of the evaluated series;
 *   funvallohilo are the third lowest doubles of the evaluated series;
 *   funvalhilolo are the second lowest doubles of the evaluated series;
 *   funvallololo are the lowest doubles of the evaluated series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   rhshihihi are the highest doubles of the linearized right hand sides;
 *   rhslohihi are the second highest doubles of the linearized rhs;
 *   rhshilohi are the third highest doubles of the linearized rhs;
 *   rhslolohi are the fourth highest doubles of the linearized rhs;
 *   rhshihilo are the fourth lowest doubles of the linearized rhs;
 *   rhslohilo are the third lowest doubles of the linearized rhs;
 *   rhshilolo are the second lowest doubles of the linearized rhs;
 *   rhslololo are the lowest doubles of the linearized right hand sides. */

void cmplx8_define_rhs
 ( int dim, int degp1,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
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
 *   Defines the right hand side for the test with the column representation
 *   of the cyclic n-roots system, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   mbrehihihi are the highest doubles of the real parts of rhs;
 *   mbrelohihi are the 2nd highest doubles of the real parts of rhs;
 *   mbrehilohi are the 3rd highest doubles of the real parts of rhs;
 *   mbrelolohi are the 4th highest doubles of the real parts of rhs;
 *   mbrehihilo are the 4th lowest doubles of the real parts of rhs;
 *   mbrelohilo are the 3rd lowest doubles of the real parts of rhs;
 *   mbrehilolo are the 2nd lowest doubles of the real parts of rhs;
 *   mbrelololo are the lowest doubles of the real parts of rhs;
 *   mbimhihihi are the highest doubles of the imag parts of rhs;
 *   mbimlohihi are the 2nd highest doubles of the imag parts of rhs;
 *   mbimhilohi are the 3rd highest doubles of the imag parts of rhs;
 *   mbimlolohi are the 4th highest doubles of the imag parts of rhs;
 *   mbimhihilo are the 4th lowest doubles of the imag parts of rhs;
 *   mbimlohilo are the 3rd lowest doubles of the imag parts of rhs;
 *   mbimhilolo are the 2nd lowest doubles of the imag parts of rhs;
 *   mbimlololo are the lowest doubles of the imag parts of rhs;
 *   funvalrehihihi are highest doubles of real parts of the evaluations;
 *   funvalrelohihi are 2nd highest doubles of real parts of evaluations;
 *   funvalrehilohi are 3rd highest doubles of real parts of evaluations;
 *   funvalrelolohi are 4th highest doubles of real parts of evaluations;
 *   funvalrehihilo are 4th lowest doubles of real parts of evaluations;
 *   funvalrelohilo are 3rd lowest doubles of real parts of evaluations;
 *   funvalrehilolo are 2nd lowest doubles of real parts of evaluations;
 *   funvalrelololo are lowest doubles of real parts of evaluations;
 *   funvalimhihihi are highest doubles of imag parts of evaluations;
 *   funvalimlohihi are 2nd highest doubles of imag parts of evaluations;
 *   funvalimhilohi are 3rd highest doubles of imag parts of evaluations;
 *   funvalimlolohi are 4th highest doubles of imag parts of evaluations;
 *   funvalimhihilo are 4th lowest doubles of imag parts of evaluations;
 *   funvalimlohilo are 3rd lowest doubles of imag parts of evaluations;
 *   funvalimhilolo are 2nd lowest doubles of imag parts of evaluations;
 *   funvalimlololo are lowest doubles of imag parts of evaluations;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   rhsrehihihi are highest doubles of real parts of the linearized rhs;
 *   rhsrelohihi are 2nd highest doubles of real parts of the linearized rhs;
 *   rhsrehilohi are 3rd highest doubles of real parts of the linearized rhs;
 *   rhsrelolohi are 4th highest doubles of real parts of the linearized rhs;
 *   rhsrehihilo are 4th lowest doubles of real parts of the linearized rhs;
 *   rhsrelohilo are 3rd lowest doubles of real parts of the linearized rhs;
 *   rhsrehilolo are 2nd lowest doubles of real parts of the linearized rhs;
 *   rhsrelololo are lowest doubles of real parts of the linearized rhs;
 *   rhsimhihihi are highest doubles of imag parts of the linearized rhs;
 *   rhsimlohihi are 2nd highest doubles of imag parts of the linearized rhs;
 *   rhsimhilohi are 3rd highest doubles of imag parts of the linearized rhs;
 *   rhsimlolohi are 4th highest doubles of imag parts of the linearized rhs;
 *   rhsimhihilo are 4th lowest doubles of imag parts of the linearized rhs;
 *   rhsimlohilo are 3rd lowest doubles of imag parts of the linearized rhs;
 *   rhsimhilolo are 2nd lowest doubles of imag parts of the linearized rhs;
 *   rhsimlololo are lowest doubles of imag parts of the linearized rhs. */

void dbl8_map_evaldiff_output
 ( int dim, int deg,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Maps the output of the evaluation and differentation to the
 *   real function values and the values of the Jacobian matrix.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   output    output[i][j][k] is the k-th coefficient of the series
 *             of the derivative of the i-th polynomial with respect
 *             to the variable j;
 *   funvalhihihi has space for the dim series of the function values;
 *   funvallohihi has space for the dim series of the function values;
 *   funvalhilohi has space for the dim series of the function values;
 *   funvallolohi has space for the dim series of the function values;
 *   funvalhihilo has space for the dim series of the function values;
 *   funvallohilo has space for the dim series of the function values;
 *   funvalhilolo has space for the dim series of the function values;
 *   funvallololo has space for the dim series of the function values;
 *   jacvalhihihi has space for the matrix series of the Jacobian.
 *   jacvallohihi has space for the matrix series of the Jacobian.
 *   jacvalhilohi has space for the matrix series of the Jacobian.
 *   jacvallolohi has space for the matrix series of the Jacobian.
 *   jacvalhihilo has space for the matrix series of the Jacobian.
 *   jacvallohilo has space for the matrix series of the Jacobian.
 *   jacvalhilolo has space for the matrix series of the Jacobian.
 *   jacvallololo has space for the matrix series of the Jacobian.
 *
 * ON RETURN :
 *   funvalhihihi has the highest doubles the function values;
 *   funvallohihi has the 2nd highest doubles the function values;
 *   funvalhilohi has the 3rd highest doubles the function values;
 *   funvallolohi has the 4th highest doubles the function values;
 *   funvalhihilo has the 4th lowest doubles the function values;
 *   funvallohilo has the 3rd lowest doubles the function values;
 *   funvalhilolo has the 2nd lowest doubles the function values;
 *   funvallololo has the lowest doubles the function values;
 *   jacvalhihihi has the highest doubles of the Jacobian matrix series;
 *   jacvallohihi has the 2nd highest doubles of the Jacobian matrix series;
 *   jacvalhilohi has the 3rd highest doubles of the Jacobian matrix series;
 *   jacvallolohi has the 4th highest doubles of the Jacobian matrix series;
 *   jacvalhihilo has the 4th lowest doubles of the Jacobian matrix series;
 *   jacvallohilo has the 3rd lowest doubles of the Jacobian matrix series;
 *   jacvalhilolo has the 2nd lowest doubles of the Jacobian matrix series;
 *   jacvallololo has the lowest doubles of the Jacobian matrix series. */

void cmplx8_map_evaldiff_output
 ( int dim, int deg,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Maps the output of the evaluation and differentation to the
 *   complex function values and the values of the Jacobian matrix.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   outputrehihihi has the highest doubles of the real parts of output;
 *   outputrelohihi has the 2nd highest doubles of the real parts of output;
 *   outputrehilohi has the 3rd highest doubles of the real parts of output;
 *   outputrelolohi has the 4th highest doubles of the real parts of output;
 *   outputrehihilo has the 4th lowest doubles of the real parts of output;
 *   outputrelohilo has the 3rd lowest doubles of the real parts of output;
 *   outputrehilolo has the 2nd lowest doubles of the real parts of output;
 *   outputrelololo has the lowest doubles of the real parts of output;
 *   funvalrehihihi has space for the dim series of the function values;
 *   funvalrelohihi has space for the dim series of the function values;
 *   funvalrehilohi has space for the dim series of the function values;
 *   funvalrelolohi has space for the dim series of the function values;
 *   funvalrehihilo has space for the dim series of the function values;
 *   funvalrelohilo has space for the dim series of the function values;
 *   funvalrehilolo has space for the dim series of the function values;
 *   funvalrelololo has space for the dim series of the function values;
 *   funvalimhihihi has space for the dim series of the function values;
 *   funvalimlohihi has space for the dim series of the function values;
 *   funvalimhilohi has space for the dim series of the function values;
 *   funvalimlolohi has space for the dim series of the function values;
 *   funvalimhihilo has space for the dim series of the function values;
 *   funvalimlohilo has space for the dim series of the function values;
 *   funvalimhilolo has space for the dim series of the function values;
 *   funvalimlololo has space for the dim series of the function values;
 *   jacvalrehihihi has space for the matrix series of the Jacobian;
 *   jacvalrelohihi has space for the matrix series of the Jacobian;
 *   jacvalrehilohi has space for the matrix series of the Jacobian;
 *   jacvalrelolohi has space for the matrix series of the Jacobian;
 *   jacvalrehihilo has space for the matrix series of the Jacobian;
 *   jacvalrelohilo has space for the matrix series of the Jacobian;
 *   jacvalrehilolo has space for the matrix series of the Jacobian;
 *   jacvalrelololo has space for the matrix series of the Jacobian;
 *   jacvalimhihihi has space for the matrix series of the Jacobian;
 *   jacvalimlohihi has space for the matrix series of the Jacobian;
 *   jacvalimhilohi has space for the matrix series of the Jacobian;
 *   jacvalimlolohi has space for the matrix series of the Jacobian;
 *   jacvalimhihilo has space for the matrix series of the Jacobian;
 *   jacvalimlohilo has space for the matrix series of the Jacobian;
 *   jacvalimhilolo has space for the matrix series of the Jacobian;
 *   jacvalimlololo has space for the matrix series of the Jacobian;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihihi has the highest doubles of the real parts
 *             of the function values;
 *   funvalrelohihi has the second highest doubles of the real parts
 *             of the function values;
 *   funvalrehilohi has the third highest doubles of the real parts
 *             of the function values;
 *   funvalrelolohi has the fourth highest doubles of the real parts
 *             of the function values;
 *   funvalrehihilo has the fourth lowest doubles of the real parts
 *             of the function values;
 *   funvalrelohilo has the third lowest doubles of the real parts
 *             of the function values;
 *   funvalrehilolo has the second lowest doubles of the real parts
 *             of the function values;
 *   funvalrelololo has the lowest doubles of the real parts
 *             of the function values;
 *   funvalimhihihi has the highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlohihi has the second highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimhilohi has the third highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlolohi has the fourth highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimhihilo has the fourth lowest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlohilo has the third lowest doubles of the imaginary parts
 *             of the function values;
 *   funvalimhilolo has the second lowest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlololo has the lowest doubles of the imaginary parts
 *             of the function values;
 *   jacvalrehihihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelohihi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehilohi_h are the third highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelolohi_h are the fourth highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehihilo_h are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelohilo_h are the third lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehilolo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelololo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhihihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlohihi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhilohi_h are the third highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlolohi_h are the fourth highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhihilo_h are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlohilo_h are the third lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhilolo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlololo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehihihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelohihi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehilohi_d are the third highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelolohi_d are the fourth highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehihilo_d are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelohilo_d are the third lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehilolo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelololo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhihihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlohihi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhilohi_d are the third highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlolohi_d are the fourth highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhihilo_d are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlohilo_d are the third lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhilolo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlololo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device. */

#endif
