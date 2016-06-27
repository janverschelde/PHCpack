with Generic_Polynomials;
with DoblDobl_Dense_Series_Ring;

package DoblDobl_Series_Polynomials is 
  new Generic_Polynomials(DoblDobl_Dense_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with double double complex coefficients.
