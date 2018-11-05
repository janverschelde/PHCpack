with Standard_Complex_Series_Ring;
with Standard_CSeries_Polynomials;
with Generic_Polynomial_Systems;

package Standard_CSeries_Poly_Systems is
  new Generic_Polynomial_Systems(Standard_Complex_Series_Ring,
                                 Standard_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines systems of polynomials in several variables with coefficients
--   as series of double precision complex numbers.
