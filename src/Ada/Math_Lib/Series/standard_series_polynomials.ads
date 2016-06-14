with Generic_Polynomials;
with Standard_Dense_Series_Ring;

package Standard_Series_Polynomials is 
  new Generic_Polynomials(Standard_Dense_Series_Ring);

-- DESCRIPTION :
--   Defines the polynomials in several variables, where the coefficients
--   are truncated power series with double precision complex coefficients.
