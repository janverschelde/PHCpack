// The file dbl4_systems_host.h specifies functions to evaluate and
// differentiate monomials with common factors in quad double precision.

#ifndef __dbl4_systems_host_h__
#define __dbl4_systems_host_h__

void CPU_dbl4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series, on real data.
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
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   2nd highest doubles of the coefficients of the monomials;
 *   cffhilo   2nd lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   acchihi   space to accumulate one power series of degree deg;
 *   acclohi   space to accumulate one power series of degree deg;
 *   acchilo   space to accumulate one power series of degree deg;
 *   acclolo   space to accumulate one power series of degree deg;
 *   inputhihi are the highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohi are the second highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilo are the second lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlolo are the lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   outputhihi has space for the highest doubles of the output;
 *   outputlohi has space for the second highest doubles of the output;
 *   outputhilo has space for the second lowest doubles of the output;
 *   outputlolo has space for the lowest doubles of the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffhihi   highest doubles of the evaluated common factors;
 *   cfflohi   second highest doubles of the evaluated common factors;
 *   cffhilo   second lowest doubles of the evaluated common factors;
 *   cfflolo   lowest doubles of the evaluated common factors;
 *   outputhihi are the highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohi are the second highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilo are the second lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlolo are the lowest doubles of the evaluated
 *             and differentiated monomials,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void CPU_cmplx4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi, double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi, double *accimhilo, double *accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi, 
   double **inputimhilo, double **inputimlolo, 
   double ***outputrehihi, double ***outputrelohi, 
   double ***outputrehilo, double ***outputrelolo, 
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo, int vrblvl );
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
 *   cffrehihi are the highest doubles of real parts of coefficients;
 *   cffrelohi are the second highest doubles of real parts of coefficients;
 *   cffrehilo are the second lowest doubles of real parts of coefficients;
 *   cffrelolo are the lowest doubles of real parts of coefficients;
 *   cffimhihi are the highest doubles of imaginary parts of coefficients;
 *   cffimlohi are the second highest doubles of imag parts of coefficients;
 *   cffimhilo are the second lowest doubles of imag parts of coefficients;
 *   cffimlolo are the lowest doubles of imaginary parts of coefficients;
 *   accrehihi has space to accumulate one power series of degree deg;
 *   accrelohi has space to accumulate one power series of degree deg;
 *   accrehilo has space to accumulate one power series of degree deg;
 *   accrelolo has space to accumulate one power series of degree deg;
 *   accimhihi has space to accumulate one power series of degree deg;
 *   accimlohi has space to accumulate one power series of degree deg;
 *   accimhilo has space to accumulate one power series of degree deg;
 *   accimlolo has space to accumulate one power series of degree deg;
 *   inputrehihi are the highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelohi are the second highest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrehilo are the second lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputrelolo are the lowest doubles of the real parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhihi are the highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlohi are the second highest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimhilo are the second lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   inputimlolo are the lowest doubles of the imaginary parts of
 *             coefficients of the series of degree deg, for dim variables;
 *   outputrehihi has space for the output;
 *   outputrelohi has space for the output;
 *   outputrehilo has space for the output;
 *   outputrelolo has space for the output;
 *   outputimhihi has space for the output;
 *   outputimlohi has space for the output;
 *   outputimhilo has space for the output;
 *   outputimlolo has space for the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffrehihi are the highest doubles of real parts of the common factors;
 *   cffrelohi are the 2nd highest doubles of real parts of the common factors;
 *   cffrehilo are the 2nd lowest doubles of real parts of the common factors;
 *   cffrelolo are the lowest doubles of real parts of the common factors;
 *   cffimhihi are the highest doubles of imag parts of the common factors;
 *   cffimlohi are the 2nd highest doubles of imag parts of the common factors;
 *   cffimhilo are the 2nd lowest doubles of imag parts of the common factors;
 *   cffimlolo are the lowest doubles of imag parts of the common factors;
 *   outputrehihi are the highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehihi[i][dim] is the
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehihi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelohi are the 2nd highest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelohi[i][dim] is the 2nd 
 *             highest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrehilo are the 2nd lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehilo[i][dim] is the 2nd
 *             lowest double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelolo are the lowest doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelo[i][dim] is the lowest
 *             double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelolo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhihi are the highest doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimhihi[i][dim] is the
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhihi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlohi are the 2nd highest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimlohi[i][dim] is the 2nd
 *             highest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlohi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhilo are the 2nd lowest doubles of the imag parts of evaluated
 *             and differentiated monomials, outputimhilo[i][dim] is the 2nd
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhilo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlolo are the lowest doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimlolo[i][dim] is the
 *             lowest double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlolo[i][idx[k]] is the derivative w.r.t. idx[k]. */

void CPU_dbl4_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double **acchihi, double **acclohi, double **acchilo, double **acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl );
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
 *   cffhihi   highest double coefficients of the monomials, cffhihi[i][j]
 *             are the coefficients of the j-th monomial in the i-th column;
 *   cfflohi   2nd highest double coefficients of the monomials, cfflohi[i][j]
 *             are the coefficients of the j-th monomial in the i-th column;
 *   cffhilo   2nd lowest double coefficients of the monomials, cffhilo[i][j]
 *             are the coefficients of the j-th monomial in the i-th column;
 *   cfflolo   lowest double coefficients of the monomials, cfflolo[i][j]
 *             are the coefficients of the j-th monomial in the i-th column;
 *   acchihi   space to accumulate dim+1 power series of degree deg;
 *   acclohi   space to accumulate dim+1 power series of degree deg;
 *   acchilo   space to accumulate dim+1 power series of degree deg;
 *   acclolo   space to accumulate dim+1 power series of degree deg;
 *   inputhihi are the highest double coefficients of the input;
 *   inputlohi are the second highest double coefficients of the input;
 *   inputhilo are the second lowest double coefficients of the input;
 *   inputlolo are the lowest double coefficients of the input;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihi are the highest doubles of evaluations;
 *   funvallohi are the second highest doubles of evaluations;
 *   funvalhilo are the second lowest doubles of evaluations;
 *   funvallolo are the lowest doubles of evaluations;
 *   jacvalhihi are the highest doubles of all derivatives;
 *   jacvallohi are the second highest doubles of all derivatives;
 *   jacvalhilo are the second lowest doubles of all derivatives;
 *   jacvallolo are the lowest doubles of all derivatives. */

void CPU_cmplx4_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **accrehihi, double **accrelohi,
   double **accrehilo, double **accrelolo,
   double **accimhihi, double **accimlohi,
   double **accimhilo, double **accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl );
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
 *   cffrehihi are the highest double real parts of the coefficients,
 *             cffrehihi[i][j] of the j-th mononmial in the i-th column;
 *   cffrelohi are the 2nd highest double real parts of the coefficients,
 *             cffrelohi[i][j] of the j-th mononmial in the i-th column;
 *   cffrehilo are the 2nd lowest double real parts of the coefficients,
 *             cffrehilo[i][j] of the j-th mononmial in the i-th column;
 *   cffrelolo are the lowest double real parts of the coefficients,
 *             cffrelolo[i][j] of the j-th mononmial in the i-th column;
 *   cffimhihi are the highest double imaginary parts of the coefficients,
 *             cffimhihi[i][j] of the j-th mononmial in the i-th column;
 *   cffimlohi are the 2nd highest double imaginary parts of the coefficients,
 *             cffimlohi[i][j] of the j-th mononmial in the i-th column;
 *   cffimhilo are the 2nd lowest double imaginary parts of the coefficients,
 *             cffimhilo[i][j] of the j-th mononmial in the i-th column;
 *   cffimlolo are the lowest double imaginary parts of the coefficients,
 *             cffimlolo[i][j] of the j-th mononmial in the i-th column;
 *   accrehihi has space to accumulate dim+1 power series of degree deg;
 *   accrelohi has space to accumulate dim+1 power series of degree deg;
 *   accrelohi has space to accumulate dim+1 power series of degree deg;
 *   accrelolo has space to accumulate dim+1 power series of degree deg;
 *   accimhihi has space to accumulate dim+1 power series of degree deg;
 *   accimlohi has space to accumulate dim+1 power series of degree deg;
 *   accimlohi has space to accumulate dim+1 power series of degree deg;
 *   accimlolo has space to accumulate dim+1 power series of degree deg;
 *   inputrehihi are the highest double real parts of the input;
 *   inputrelohi are the second highest double real parts of the input;
 *   inputrehilo are the second lowest double real parts of the input;
 *   inputrelolo are the lowest double real parts of the input;
 *   inputimhihi are the highest double imaginary parts of the input;
 *   inputimlohi are the second highest double imaginary parts of the input;
 *   inputimhilo are the second lowest double imaginary parts of the input;
 *   inputimlolo are the lowest double imaginary parts of the input;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihi are highest doubles of real parts of the evaluations;
 *   funvalrelohi are 2nd highest doubles of real parts of the evaluations;
 *   funvalrehilo are 2nd lowest doubles of real parts of the evaluations;
 *   funvalrelolo are lowest doubles of real parts of the evaluations;
 *   funvalimhihi are highest doubles of imag parts of the evaluations;
 *   funvalimlohi are 2nd highest doubles of imag parts of the evaluations;
 *   funvalimhilo are 2nd lowest doubles of imag parts of the evaluations;
 *   funvalimlolo are lowest doubles of imag parts of the evaluations;
 *   jacvalrehihi are highest doubles of real parts of all derivatives;
 *   jacvalrelohi are 2nd highest doubles of real parts of all derivatives;
 *   jacvalrehilo are 2nd lowest doubles of real parts of all derivatives;
 *   jacvalrelolo are lowest doubles of real parts of all derivatives;
 *   jacvalimhihi are highest doubles of imag parts of all derivatives;
 *   jacvalimlohi are 2nd highest doubles of imag parts of all derivatives;
 *   jacvalimhilo are 2nd lowest doubles of imag parts of all derivatives;
 *   jacvalimlolo are lowest doubles of imag parts of all derivatives. */

void dbl4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double damper,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi, 
   double **funvalhilo, double **funvallolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl );
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
 *   mbhihi    highest doubles of the right hand side of monomial system;
 *   mblohi    second highest doubles of the right hand side;
 *   mbhilo    second lowest doubles of the right hand side;
 *   mblolo    lowest doubles of the right hand side;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputhihi are the highest doubles of the output
 *             of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputlohi are the second highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilo are the second lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlolo are the lowest doubles
 *             of the output of the evaluation and differentiation;
 *   funvalhihi has space allocated for dim power series;
 *   funvallohi has space allocated for dim power series;
 *   funvalhilo has space allocated for dim power series;
 *   funvallolo has space allocated for dim power series;
 *   rhshihi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslolo   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalhihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihi has the highest doubles of output[i][dim];
 *   funvallohi has the second highest doubles of output[i][dim];
 *   funvallohi has the second lowest doubles of output[i][dim];
 *   funvallolo has the lowest doubles of output[i][dim];
 *   rhshihi   the linearized right hand side are the highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohi   the linearized right hand side are the 2nd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilo   the linearized right hand side are the 2nd lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslolo   the linearized right hand side are the lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   jacvalhihi is a series with matrices as coefficients, of highest doubles,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi is a series with matrices as coefficients, of second highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilo is a series with matrices as coefficients, of second lowest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallolo is a series with matrices as coefficients, of lowest doubles,
 *             the leading coefficient is the Jacobian matrix. */

void cmplx4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double damper,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl );
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
 *   mbrehihi  highest real parts of the right hand side series;
 *   mbrelohi  second highest real parts of the right hand side series;
 *   mbrehilo  second lowest real parts of the right hand side series;
 *   mbrelolo  lowest real parts of the right hand side series;
 *   mbimhihi  highest imaginary parts of the right hand side series;
 *   mbimlohi  second highest imaginary parts of the right hand side series;
 *   mbimhilo  second lowest imaginary parts of the right hand side series;
 *   mbimlolo  lowest imaginary parts of the right hand side series;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputrehihi are the highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehihi[i][dim] is the power series value, 
 *             outputrehihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelohi are the 2nd highest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelohi[i][dim] is the power series value, 
 *             outputrelohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrehilo are the 2nd lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehilo[i][dim] is the power series value, 
 *             outputrehilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelolo are the lowest doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelolo[i][dim] is the power series value, 
 *             outputrelolo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhihi are the highest doubles of the imaginary parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhihi[i][dim] is the power series value, 
 *             outputimhihi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlohi are the 2nd highest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlohi[i][dim] is the power series value, 
 *             outputimlohi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhilo are the 2nd lowest doubles of the imag parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhilo[i][dim] is the power series value, 
 *             outputimhilo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlolo are the lowest doubles of the imaginary parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlolo[i][dim] is the power series value, 
 *             outputimlolo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   funvalrehihi has space allocated for dim power series;
 *   funvalrelohi has space allocated for dim power series;
 *   funvalrehilo has space allocated for dim power series;
 *   funvalrelolo has space allocated for dim power series;
 *   funvalimhihi has space allocated for dim power series;
 *   funvalimlohi has space allocated for dim power series;
 *   funvalimhilo has space allocated for dim power series;
 *   funvalimlolo has space allocated for dim power series;
 *   rhsrehihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrehilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelolo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlohi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlolo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalrehihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrehilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihi are the highest doubles of the real parts in output[i][dim];
 *   funvalrelohi are the second highest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrehilo are the second lowest doubles of the real parts
 *             of the output[i][dim];
 *   funvalrelolo are the lowest doubles of the real parts in output[i][dim];
 *   funvalimhihi are the highest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimlohi are the second highest doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimhilo are the second lowest doubles of the imaginary parts
 *             of the output[i][dim];
 *   funvalimlolo are the lowest doubles of the imaginary parts
 *             of the output[i][dim];
 *   rhsrehihi are the highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelohi are the second highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrehilo are the second lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsrelolo are the lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t;
 *   rhsimhihi are the highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlohi are the second highest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimhilo are the second lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   rhsimlolo are the lowest doubles of imaginary parts of
 *             the linearized right hand side, subtracted by 1 and added by t;
 *   jacvalrehihi are the highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelohi are the second highest doubles of the real parts
 *             of the matrix series;
 *   jacvalrehilo are the second lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalrelolo are the lowest doubles of the real parts
 *             of the matrix series;
 *   jacvalimhihi are the highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimlohi are the second highest doubles of imaginary parts
 *             of the matrix series;
 *   jacvalimhilo are the second lowest doubles imaginary parts
 *             of the matrix series;
 *   jacvalimlolo are the lowest doubles imaginary parts
 *             of the matrix series. */

void dbl4_define_rhs
 ( int dim, int degp1,
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double **rhshihi, double **rhslohi,
   double **rhshilo, double **rhslolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the right hand side for the test with the column representation
 *   of the cyclic n-roots system, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   mbhihi    highest doubles of the right hand side series;
 *   mblohi    second highest doubles of the right hand side series;
 *   mbhilo    second lowest doubles of the right hand side series;
 *   mblolo    lowest doubles of the right hand side series;
 *   funvalhihi are the highest doubles of the evaluated series;
 *   funvallohi are the second highest doubles of the evaluated series;
 *   funvalhilo are the second lowest doubles of the evaluated series;
 *   funvallolo are the lowest doubles of the evaluated series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   rhshihi   highest doubles of the linearized right hand sides;
 *   rhslohi   second highest doubles of the linearized right hand sides;
 *   rhshilo   second lowest doubles of the linearized right hand sides;
 *   rhslolo   lowest doubles of the linearized right hand sides. */

void cmplx4_define_rhs
 ( int dim, int degp1,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the right hand side for the test with the column representation
 *   of the cyclic n-roots system, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   mbrehihi  highest doubles of the real parts of the right hand side;
 *   mbrelohi  2nd highest doubles of the real parts of the right hand side;
 *   mbrehilo  2nd lowest doubles of the real parts of the right hand side;
 *   mbrelolo  lowest doubles of the real parts of the right hand side;
 *   mbimhihi  highest doubles of the imag parts of the right hand side;
 *   mbimlohi  2nd highest doubles of the imag parts of the right hand side;
 *   mbimhilo  2nd lowest doubles of the imag parts of the right hand side;
 *   mbimlolo  lowest doubles of the imag parts of the right hand side;
 *   funvalrehihi are highest doubles of the real parts of the evaluations;
 *   funvalrelohi are 2nd highest doubles of the real parts of the evaluations;
 *   funvalrehilo are 2nd lowest doubles of the real parts of the evaluations;
 *   funvalrelolo are lowest doubles of the real parts of the evaluations;
 *   funvalimhihi are highest doubles of the imag parts of the evaluations;
 *   funvalimlohi are 2nd highest doubles of the imag parts of the evaluations;
 *   funvalimhilo are 2nd lowest doubles of the imag parts of the evaluations;
 *   funvalimlolo are lowest doubles of the imag parts of the evaluations;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   rhsrehihi are highest doubles of real parts of the linearized rhs;
 *   rhsrelohi are 2nd highest doubles of real parts of the linearized rhs;
 *   rhsrehilo are 2nd lowest doubles of real parts of the linearized rhs;
 *   rhsrelolo are lowest doubles of real parts of the linearized rhs;
 *   rhsimhihi are highest doubles of imag parts of the linearized rhs;
 *   rhsimlohi are 2nd highest doubles of imag parts of the linearized rhs;
 *   rhsimhilo are 2nd lowest doubles of imag parts of the linearized rhs;
 *   rhsimlolo are lowest doubles of imag parts of the linearized rhs. */

void dbl4_map_evaldiff_output
 ( int dim, int deg,
    double ***outputhihi, double ***outputlohi,
    double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl );
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
 *   funvalhihi has space for the dim series of the function values;
 *   funvallohi has space for the dim series of the function values;
 *   funvalhilo has space for the dim series of the function values;
 *   funvallolo has space for the dim series of the function values;
 *   jacvalhihi has space for the matrix series of the Jacobian.
 *   jacvallohi has space for the matrix series of the Jacobian.
 *   jacvalhilo has space for the matrix series of the Jacobian.
 *   jacvallolo has space for the matrix series of the Jacobian.
 *
 * ON RETURN :
 *   funvalhihi has the highest doubles the function values;
 *   funvallohi has the 2nd highest doubles the function values;
 *   funvalhilo has the 2nd lowest doubles the function values;
 *   funvallolo has the lowest doubles the function values;
 *   jacvalhihi has the highest doubles of the Jacobian matrix series;
 *   jacvallohi has the 2nd highest doubles of the Jacobian matrix series;
 *   jacvalhilo has the 2nd lowest doubles of the Jacobian matrix series;
 *   jacvallolo has the lowest doubles of the Jacobian matrix series. */

void cmplx4_map_evaldiff_output
 ( int dim, int deg,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl );
/*
 * DESCRIPTION :
 *   Maps the output of the evaluation and differentation to the
 *   complex function values and the values of the Jacobian matrix.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   outputrehihi has the highest doubles of the real parts of the output;
 *   outputrelohi has the 2nd highest doubles of the real parts of the output;
 *   outputrehilo has the 2nd lowest doubles of the real parts of the output;
 *   outputrelolo has the lowest doubles of the real parts of the output;
 *   funvalrehihi has space for the dim series of the function values;
 *   funvalrelohi has space for the dim series of the function values;
 *   funvalrehilo has space for the dim series of the function values;
 *   funvalrelolo has space for the dim series of the function values;
 *   funvalimhihi has space for the dim series of the function values;
 *   funvalimlohi has space for the dim series of the function values;
 *   funvalimhilo has space for the dim series of the function values;
 *   funvalimlolo has space for the dim series of the function values;
 *   jacvalrehihi has space for the matrix series of the Jacobian;
 *   jacvalrelohi has space for the matrix series of the Jacobian;
 *   jacvalrehilo has space for the matrix series of the Jacobian;
 *   jacvalrelolo has space for the matrix series of the Jacobian;
 *   jacvalimhihi has space for the matrix series of the Jacobian;
 *   jacvalimlohi has space for the matrix series of the Jacobian;
 *   jacvalimhilo has space for the matrix series of the Jacobian;
 *   jacvalimlolo has space for the matrix series of the Jacobian;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehihi has the highest doubles of the real parts
 *             of the function values;
 *   funvalrelohi has the second highest doubles of the real parts
 *             of the function values;
 *   funvalrehilo has the second lowest doubles of the real parts
 *             of the function values;
 *   funvalrelolo has the lowest doubles of the real parts
 *             of the function values;
 *   funvalimhihi has the highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlohi has the second highest doubles of the imaginary parts
 *             of the function values;
 *   funvalimhilo has the second lowest doubles of the imaginary parts
 *             of the function values;
 *   funvalimlolo has the lowest doubles of the imaginary parts
 *             of the function values;
 *   jacvalrehihi has the highest doubles of the real parts
 *             of the Jacobian matrix series;
 *   jacvalrelohi has the second highest doubles of the real parts
 *             of the Jacobian matrix series;
 *   jacvalrehilo has the second lowest doubles of the real parts
 *             of the Jacobian matrix series;
 *   jacvalrelolo has the lowest doubles of the real parts
 *             of the Jacobian matrix series;
 *   jacvalimhihi has the highest doubles of the imaginary parts
 *             of the Jacobian matrix series;
 *   jacvalimlohi has the second highest doubles of the imaginary parts
 *             of the Jacobian matrix series;
 *   jacvalimhilo has the second lowest doubles of the imaginary parts
 *             of the Jacobian matrix series;
 *   jacvalimlolo has the lowest doubles of the imaginary parts
 *             of the Jacobian matrix series. */

#endif
