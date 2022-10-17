// The file dbl8_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in octo double precision.

#ifndef __dbl8_monomial_systems_h__
#define __dbl8_monomial_systems_h__

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

#endif
