with Generic_Polynomials;
with DecaDobl_Complex_Series_Ring;

package DecaDobl_CSeries_Polynomials is 
  new Generic_Polynomials(DecaDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with deca double complex coefficients.
