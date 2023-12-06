with Generic_Polynomials;
with HexaDobl_Complex_Series_Ring;

package HexaDobl_CSeries_Polynomials is 
  new Generic_Polynomials(HexaDobl_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with hexa double complex coefficients.
