with Generic_Polynomials;
with OctoDobl_Complex_Series_Ring;

package OctoDobl_CSeries_Polynomials is 
  new Generic_Polynomials(OctoDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with octo double complex coefficients.
