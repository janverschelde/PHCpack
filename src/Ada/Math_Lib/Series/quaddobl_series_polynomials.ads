with Generic_Polynomials;
with QuadDobl_Dense_Series_Ring;

package QuadDobl_Series_Polynomials is 
  new Generic_Polynomials(QuadDobl_Dense_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with quad double complex coefficients.
