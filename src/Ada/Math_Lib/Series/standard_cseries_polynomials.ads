with Generic_Polynomials;
with Standard_Complex_Series_Ring;

package Standard_CSeries_Polynomials is 
  new Generic_Polynomials(Standard_Complex_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with double precision complex coefficients.
