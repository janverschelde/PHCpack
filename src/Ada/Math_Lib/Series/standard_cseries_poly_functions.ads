with Standard_Complex_Series_Ring;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;
with Generic_Polynomial_Functions;

package Standard_CSeries_Poly_Functions is
  new Generic_Polynomial_Functions(Standard_Complex_Series_Ring,
                                   Standard_Complex_Series_Vectors,
                                   Standard_CSeries_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions to evaluate a truncate power series
--   with double precision complex coefficients.
