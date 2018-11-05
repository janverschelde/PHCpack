with Generic_Polynomials;
with DoblDobl_Complex_Series_Ring;

package DoblDobl_CSeries_Polynomials is 
  new Generic_Polynomials(DoblDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with double double complex coefficients.
