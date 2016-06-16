with Standard_Dense_Series_Ring;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Generic_Polynomial_Functions;

package Standard_Series_Poly_Functions is
  new Generic_Polynomial_Functions(Standard_Dense_Series_Ring,
                                   Standard_Dense_Series_Vectors,
                                   Standard_Series_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with double precision complex coefficients.
