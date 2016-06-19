with Standard_Dense_Series_Ring;
with Standard_Series_Polynomials;
with Generic_Polynomial_Systems;

package Standard_Series_Poly_Systems is
  new Generic_Polynomial_Systems(Standard_Dense_Series_Ring,
                                 Standard_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with polynomials with coefficients
--   as series of complex numbers.
