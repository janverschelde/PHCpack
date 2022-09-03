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

void dbl8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
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

#endif
