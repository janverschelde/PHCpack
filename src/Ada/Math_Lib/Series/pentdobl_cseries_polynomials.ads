with Generic_Polynomials;
with PentDobl_Complex_Series_Ring;

package PentDobl_CSeries_Polynomials is 
  new Generic_Polynomials(PentDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with penta double complex coefficients.
